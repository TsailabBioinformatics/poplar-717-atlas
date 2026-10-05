#!/usr/bin/env python3
"""A static file server that honours HTTP Range, for local checks.

`python3 -m http.server` does NOT implement Range: it ignores the header and answers 200 with
the whole file. That is exactly the failure the expression loader now refuses, so testing
against it would exercise the fallback path and never the real one -- and worse, a client that
did NOT refuse would inflate the first gzip member in the blob and render another gene's
expression profile without erroring.

GitHub Pages does honour Range (verified: real 206, identity-encoded), so this is only about
making `check.mjs` test what production actually does.

    python3 scripts/serve_ranges.py 8934
"""
import functools, http.server, os, re, socketserver, sys

RANGE_RE = re.compile(r"^bytes=(\d*)-(\d*)$")


class RangeHandler(http.server.SimpleHTTPRequestHandler):
    def do_GET(self):                                     # noqa: N802
        rng = self.headers.get("Range")
        if not rng:
            return super().do_GET()
        path = self.translate_path(self.path.split("?", 1)[0])
        if not os.path.isfile(path):
            return super().do_GET()
        m = RANGE_RE.match(rng.strip())
        if not m:
            self.send_error(416, "malformed Range")
            return
        size = os.path.getsize(path)
        start_s, end_s = m.group(1), m.group(2)
        if start_s == "":                                  # suffix range: last N bytes
            length = int(end_s or 0)
            start, end = max(0, size - length), size - 1
        else:
            start = int(start_s)
            end = int(end_s) if end_s else size - 1
        if start >= size or end < start:
            self.send_response(416)
            self.send_header("Content-Range", f"bytes */{size}")
            self.end_headers()
            return
        end = min(end, size - 1)
        n = end - start + 1
        self.send_response(206)
        self.send_header("Content-Type", self.guess_type(path))
        self.send_header("Content-Range", f"bytes {start}-{end}/{size}")
        self.send_header("Content-Length", str(n))
        self.send_header("Accept-Ranges", "bytes")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        with open(path, "rb") as fh:
            fh.seek(start)
            self.wfile.write(fh.read(n))

    def end_headers(self):
        # Pages sends accept-ranges on ordinary responses too; mirror it so a client that
        # probes for capability sees the same thing locally as in production.
        if "Accept-Ranges" not in self._headers_buffer_names():
            self.send_header("Accept-Ranges", "bytes")
        super().end_headers()

    def _headers_buffer_names(self):
        return [h.split(b":")[0].decode("latin-1") for h in getattr(self, "_headers_buffer", [])
                if b":" in h]

    def log_message(self, *a):                            # quiet
        pass


if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8934
    root = sys.argv[2] if len(sys.argv) > 2 else os.getcwd()
    handler = functools.partial(RangeHandler, directory=root)
    socketserver.ThreadingTCPServer.allow_reuse_address = True
    with socketserver.ThreadingTCPServer(("127.0.0.1", port), handler) as httpd:
        print(f"serving {root} on 127.0.0.1:{port} with Range support")
        httpd.serve_forever()
