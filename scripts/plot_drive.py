#!/usr/bin/env python3
"""Publish the local CheMPAS plot outbox to an SMB share or Google Drive."""
import argparse
import configparser
import getpass
import json
import os
from pathlib import Path
import shutil
import socket
import subprocess
import tempfile
from urllib.parse import quote

REMOTE = "chempas-plots"
SMB_REMOTE = "chempas-mac"
CONFIG = Path.home() / ".config/rclone/rclone.conf"
SETTINGS = CONFIG.with_name("chempas-plots.json")
OUTBOX = Path.home() / "Data/CheMPAS-A/plot-outbox"
DESTINATION = f"{REMOTE}:CheMPAS-A/plots"
EXTENSIONS = (".png", ".pdf", ".svg")


def read_config():
    config = configparser.ConfigParser(interpolation=None)
    config.read(CONFIG)
    return config


def write_config(config):
    CONFIG.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    CONFIG.parent.chmod(0o700)
    with tempfile.NamedTemporaryFile(mode="w", dir=CONFIG.parent, delete=False) as handle:
        config.write(handle)
        temporary = Path(handle.name)
    temporary.chmod(0o600)
    temporary.replace(CONFIG)


def default_destination():
    if SETTINGS.exists():
        return json.loads(SETTINGS.read_text())["destination"]
    return DESTINATION


def save_destination(destination):
    with tempfile.NamedTemporaryFile(mode="w", dir=CONFIG.parent, delete=False) as handle:
        json.dump({"destination": destination}, handle, indent=2)
        handle.write("\n")
        temporary = Path(handle.name)
    temporary.chmod(0o600)
    temporary.replace(SETTINGS)


def executable():
    name = shutil.which("rclone")
    if not name:
        raise SystemExit("rclone is not on PATH. Install it from https://rclone.org/install/.")
    return name


def keyring(action, host, user, password=None):
    request = {"action": action, "host": host, "user": user}
    if password is not None:
        request["password"] = password
    result = subprocess.run(
        ["/usr/bin/python3", str(Path(__file__).resolve().with_name("plot_keyring.py"))],
        input=json.dumps(request), text=True, capture_output=True)
    if result.returncode:
        raise SystemExit("Could not access the desktop login keyring. Unlock it and try again.")
    return json.loads(result.stdout)


def smb_environment(password):
    secret = subprocess.run([executable(), "obscure", "-"], input=password + "\n",
                            text=True, capture_output=True, check=True).stdout.strip()
    environment = os.environ.copy()
    environment["RCLONE_SMB_PASS"] = secret
    return environment


def smb_password(host, user, share, gui):
    if not gui:
        return getpass.getpass(f"Mac login password for {user}@{host}: ")
    import tkinter as tk
    from tkinter import simpledialog

    root = tk.Tk()
    root.withdraw()
    root.attributes("-topmost", True)
    print("Password dialog opened on this computer; enter your Mac login password there.", flush=True)
    try:
        password = simpledialog.askstring(
            "CheMPAS plots — rclone SMB login",
            f"Mac login password for {user}@{host}\n"
            f"Shared folder: {share}\n\n"
            "Saved in this computer's login keyring for future plot uploads.\n"
            "The password will not be written into rclone.conf or the repository.",
            show="*", parent=root)
    finally:
        root.destroy()
    if password is None:
        raise SystemExit("Connection cancelled; configuration was not changed.")
    return password


def connect_smb(host, user, share, gui):
    if not host or not user or any(c in host + user + share for c in "\r\n"):
        raise SystemExit("Supply a hostname and Mac login username without line breaks.")
    if not share or share in (".", "..") or "/" in share or "\\" in share:
        raise SystemExit("Supply the shared folder name, such as CheMPAS-A-Plots.")
    # Resolve Bonjour names through the OS before rclone's DNS-only resolver.
    try:
        host = socket.gethostbyname(host)
    except OSError:
        raise SystemExit("Could not resolve the Mac hostname. Use its LAN IP address.") from None
    config = read_config()
    if config.has_section(SMB_REMOTE) and config[SMB_REMOTE].get("type") != "smb":
        raise SystemExit(f"Existing {SMB_REMOTE} remote is not SMB; leaving it unchanged.")
    password = smb_password(host, user, share, gui)
    if not password or "\n" in password or "\r" in password:
        raise SystemExit("A nonempty Mac login password without line breaks is required.")
    environment = smb_environment(password)
    config[SMB_REMOTE] = {"type": "smb", "host": host, "user": user,
                          "description": "Password in desktop login keyring; use chempas-plots upload"}
    destination = f"{SMB_REMOTE}:{share}"
    # Validate using a private temporary config before replacing any saved connection.
    with tempfile.TemporaryDirectory(prefix="chempas-smb-") as directory:
        candidate = Path(directory) / "rclone.conf"
        with candidate.open("w") as handle:
            config.write(handle)
        candidate.chmod(0o600)
        command = [executable(), "--config", str(candidate), "--contimeout", "8s", "--timeout", "15s",
                   "--retries", "1", "--low-level-retries", "1"]

        def query(path, *flags):
            result = subprocess.run(command + ["lsjson", path, *flags],
                                    text=True, capture_output=True, env=environment)
            if result.returncode:
                details = result.stderr.replace(password, "[hidden]").replace(
                    environment["RCLONE_SMB_PASS"], "[hidden]")
                raise SystemExit("SMB connection failed; configuration was not changed.\n" + details[-1600:])
            return result

        shares = json.loads(query(SMB_REMOTE + ":", "--dirs-only").stdout)
        print("Mac login succeeded.", flush=True)
        names = [item["Name"] for item in shares if item.get("IsDir")]
        matched = next((name for name in names if name.casefold() == share.casefold()), None)
        if matched:
            destination = f"{SMB_REMOTE}:{matched}"
        elif user in names:
            destination = f"{SMB_REMOTE}:{user}/{share}"
        else:
            raise SystemExit(f"Login succeeded, but {share} is not an available share.\n"
                             "Shared folders: " + ", ".join(names))
        result = query(destination, "--stat")
        if not json.loads(result.stdout).get("IsDir"):
            raise SystemExit("The shared destination is not a directory; configuration was not changed.")
    keyring("store", host, user, password)
    del password, environment
    write_config(config)
    save_destination(destination)
    print(f"SMB folder access verified. Default destination: {destination}")
    print("Run: chempas-plots upload")


def connect(client_json):
    config = read_config()
    if config.has_section(REMOTE) and config[REMOTE].get("type") != "drive":
        raise SystemExit(f"Existing {REMOTE} remote is not Google Drive; leaving it unchanged.")
    if client_json:
        client = json.loads(client_json.expanduser().read_text()).get("installed", {})
        if not client.get("client_id") or not client.get("client_secret"):
            raise SystemExit("Expected the downloaded JSON for a Google OAuth Desktop app client.")
        if any("\n" in client[key] or "\r" in client[key] for key in ("client_id", "client_secret")):
            raise SystemExit("Invalid OAuth client fields.")
        secret = subprocess.run([executable(), "obscure", "-"], input=client["client_secret"] + "\n",
                                text=True, capture_output=True, check=True).stdout.strip()
        if CONFIG.exists():
            backup = CONFIG.with_name("rclone.conf.before-chempas-connect")
            shutil.copyfile(CONFIG, backup)
            backup.chmod(0o600)
        config[REMOTE] = {"type": "drive", "scope": "drive.file",
                          "client_id": client["client_id"], "client_secret": secret}
        write_config(config)
    section = config[REMOTE] if config.has_section(REMOTE) else {}
    if not section.get("client_id") or not section.get("client_secret"):
        raise SystemExit("Supply your OAuth Desktop app JSON: chempas-plots connect /path/to/client.json\n"
                         "Setup: https://rclone.org/drive/#making-your-own-client-id")
    print("Complete Google sign-in in your browser. Choose My Drive (not a Shared Drive).", flush=True)
    subprocess.run([executable(), "--config", str(CONFIG), "config", "reconnect", REMOTE + ":"], check=True)
    config = read_config()
    if not config[REMOTE].get("token"):
        raise SystemExit("Google authorization did not complete; no token was saved.")
    # Verify live access; do not print account contents or tokens.
    subprocess.run([executable(), "--config", str(CONFIG), "lsjson", REMOTE + ":", "--dirs-only"],
                   stdout=subprocess.DEVNULL, check=True)
    save_destination(DESTINATION)
    print("Google Drive authorization verified. Run: chempas-plots upload")


def upload(source, destination, dry_run):
    source = source.expanduser().resolve()
    if not source.is_dir():
        raise SystemExit(f"Plot outbox does not exist: {source}")
    files = [p for p in source.rglob("*") if p.is_file() and not p.is_symlink()
             and p.suffix.lower() in EXTENSIONS]
    if not files:
        raise SystemExit(f"No PNG, PDF, or SVG plots found in {source}")
    config = read_config()
    remote = destination.split(":", 1)[0] if ":" in destination else None
    section = config[remote] if remote and config.has_section(remote) else {}
    if section.get("type") == "drive":
        if not section.get("token"):
            raise SystemExit("Google authorization is pending. Run: chempas-plots connect /path/to/client.json")
    filters = ["--ignore-case", "--include", "*.png", "--include", "*.pdf", "--include", "*.svg"]
    command = [executable(), "--config", str(CONFIG)]
    print(f"{'Previewing' if dry_run else 'Uploading'} {len(files)} plots: {source} -> {destination}", flush=True)
    # SMB does not provide content hashes. Copy the small plot outbox in full and
    # read the files back to verify bytes, including same-size replacements.
    smb = section.get("type") == "smb"
    environment = None
    if smb and remote == SMB_REMOTE and not section.get("pass"):
        password = keyring("lookup", section["host"], section["user"])["password"]
        environment = smb_environment(password)
        del password
    copy_flags = ["--ignore-times"] if smb else ["--checksum"]
    check_flags = ["--download"] if smb else []
    subprocess.run(command + ["copy", str(source), destination] + copy_flags + filters
                   + (["--dry-run"] if dry_run else []), check=True, env=environment)
    if dry_run:
        return
    subprocess.run(command + ["check", str(source), destination, "--one-way"] + check_flags + filters,
                   check=True, env=environment)
    print("Plot contents verified.")
    if smb:
        print("SMB folder: smb://" + section["host"] + "/" + quote(destination.split(":", 1)[1], safe="/"))
    elif section.get("type") == "drive":
        result = subprocess.run(command + ["lsjson", destination, "--stat"],
                                text=True, capture_output=True, check=True)
        folder = json.loads(result.stdout)
        if folder.get("ID"):
            print("Drive folder: https://drive.google.com/drive/folders/" + folder["ID"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("status", help="Show local setup and authorization state")
    auth = commands.add_parser("connect", help="Import a Desktop OAuth client and authorize Google Drive")
    auth.add_argument("client_json", type=Path, nargs="?")
    smb = commands.add_parser("connect-smb", help="Connect a Mac/shared folder and select it for future uploads")
    smb.add_argument("host", help="Mac hostname, such as Gemini.local")
    smb.add_argument("share", help="Shared folder name, such as CheMPAS-A-Plots")
    smb.add_argument("--user", default=getpass.getuser(), help="Mac login username")
    smb.add_argument("--gui", action="store_true", help="Ask for the password in a local graphical dialog")
    publish = commands.add_parser("upload", help="Copy plots and verify their contents at the destination")
    publish.add_argument("--source", type=Path, default=OUTBOX)
    publish.add_argument("--destination", default=default_destination())
    publish.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    if args.command == "status":
        config = read_config()
        destination = default_destination()
        remote = destination.split(":", 1)[0] if ":" in destination else None
        section = config[remote] if remote and config.has_section(remote) else {}
        subprocess.run([executable(), "version"], check=True)
        print(f"Config: {CONFIG}\nPlot outbox: {OUTBOX}\nDestination: {destination}")
        if section.get("type") == "smb":
            print(f"SMB account: {section['user']}@{section['host']}")
            print("Password storage: " + ("rclone config" if section.get("pass") else "desktop login keyring"))
        elif section.get("type") == "drive":
            print("OAuth client: " + ("configured" if section.get("client_id") else "needed"))
            print("Google token: " + ("saved (not checked online)" if section.get("token") else "pending"))
    elif args.command == "connect":
        connect(args.client_json)
    elif args.command == "connect-smb":
        connect_smb(args.host, args.user, args.share, args.gui)
    else:
        upload(args.source, args.destination, args.dry_run)


if __name__ == "__main__":
    os.umask(0o077)
    try:
        main()
    except subprocess.CalledProcessError as error:
        raise SystemExit(f"rclone exited with status {error.returncode}; operation did not complete.") from None
