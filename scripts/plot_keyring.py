#!/usr/bin/python3
"""Private stdin/stdout bridge to the Linux login keyring for plot_drive.py."""
import json
import sys

import gi
gi.require_version("Secret", "1")
from gi.repository import Secret


def main():
    request = json.load(sys.stdin)
    schema = Secret.Schema.new("org.chempas.PlotPublisher", Secret.SchemaFlags.NONE,
                               {name: Secret.SchemaAttributeType.STRING
                                for name in ("application", "host", "user")})
    attributes = {"application": "chempas-plots", "host": request["host"].lower(),
                  "user": request["user"]}
    service = Secret.Service.get_sync(Secret.ServiceFlags.OPEN_SESSION, None)
    collection = Secret.Collection.for_alias_sync(service, "default", Secret.CollectionFlags.NONE, None)
    if collection is None or collection.get_object_path().endswith("/session"):
        raise RuntimeError("A persistent login keyring is required.")
    action = request["action"]
    if action == "store":
        result = Secret.password_store_sync(
            schema, attributes, collection.get_object_path(),
            f"CheMPAS plots: {request['user']}@{request['host']}", request["password"], None)
        if not result:
            raise RuntimeError("The keyring did not save the password.")
        print(json.dumps({"stored": True}))
    elif action == "lookup":
        password = Secret.password_lookup_sync(schema, attributes, None)
        if password is None:
            raise RuntimeError("No plot-sharing password is saved for this account.")
        # Only the parent helper reads this pipe; never invoke lookup in a tool log.
        print(json.dumps({"password": password}))
    elif action == "clear":
        print(json.dumps({"cleared": Secret.password_clear_sync(schema, attributes, None)}))
    else:
        raise RuntimeError("Unknown keyring action.")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        # Suppress tracebacks/payloads so credentials cannot appear in diagnostics.
        print("Login keyring operation failed. Check that the desktop keyring is unlocked.", file=sys.stderr)
        sys.exit(1)
