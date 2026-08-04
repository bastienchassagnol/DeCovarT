#!/usr/bin/env bash
# Repair home zotcurate for Zotero 8+ native citation keys (Zotero 9 OK).
# Run as bastien (fix ~/.local ownership first if needed).
set -euo pipefail
HOME_BIN="${HOME}/.local/bin"
ZOTC_DIR="${HOME}/.zotc"
ZOTERO_LIBRARY_ID="${ZOTERO_LIBRARY_ID:-7172053}"
ZOTERO_DB="${ZOTERO_DB:-${HOME}/Zotero/zotero.sqlite}"

if ! command -v uv >/dev/null 2>&1; then
  echo "ERROR: uv required" >&2
  exit 1
fi

echo "==> uv tool install upstream zotcurate"
uv tool install --force "git+https://github.com/jeetsukumaran/zotcurate.git"

echo "==> Write ${HOME_BIN}/zotc wrapper (native citationKey in zotero.sqlite)"
mkdir -p "${HOME_BIN}"
# Replace uv-tool symlink with an explicit wrapper that always passes -z.
cat > "${HOME_BIN}/zotc" <<WRAP
#!/usr/bin/env bash
set -euo pipefail
# Prefer the uv-tool environment's zotcurate (stdlib); native keys live in zotero.sqlite.
ZOTERO_DB="\${ZOTERO_DB:-${ZOTERO_DB}}"
if command -v uv >/dev/null 2>&1 && uv tool run --from zotcurate zotc -h >/dev/null 2>&1; then
  exec uv tool run --from zotcurate zotc -z "\${ZOTERO_DB}" "\$@"
fi
# Fallback: module on PYTHONPATH / site-packages
exec python3 -m zotcurate.cli -z "\${ZOTERO_DB}" "\$@"
WRAP
chmod 755 "${HOME_BIN}/zotc"

echo "==> Configure ${ZOTC_DIR}"
mkdir -p "${ZOTC_DIR}"
printf '%s\n' "${ZOTERO_LIBRARY_ID}" > "${ZOTC_DIR}/library"
printf '%s\n' "user" > "${ZOTC_DIR}/library-type"
if [[ -n "${ZOTERO_API_KEY:-}" ]]; then
  printf '%s\n' "${ZOTERO_API_KEY}" > "${ZOTC_DIR}/api-key"
elif [[ ! -f "${ZOTC_DIR}/api-key" ]]; then
  echo "ERROR: set ZOTERO_API_KEY or create ${ZOTC_DIR}/api-key" >&2
  exit 1
fi
chmod 700 "${ZOTC_DIR}"
chmod 600 "${ZOTC_DIR}/api-key" "${ZOTC_DIR}/library" "${ZOTC_DIR}/library-type"

export PATH="${HOME_BIN}:${PATH}"
hash -r 2>/dev/null || true
echo "==> Smoke test (copy DB if locked by desktop Zotero)"
command -v zotc
# Avoid lock: copy then point -z at the copy for the smoke test only
TMP_DB="$(mktemp /tmp/zotero-smoke-XXXXXX.sqlite)"
cp -f "${ZOTERO_DB}" "${TMP_DB}"
zotc -z "${TMP_DB}" config || ZOTERO_DB="${TMP_DB}" zotc config
echo "keys: $(zotc -z "${TMP_DB}" keys list 2>/dev/null | wc -l | tr -d ' ')"
rm -f "${TMP_DB}"
echo "Done. Day-to-day: zotc -z ${ZOTERO_DB} … (quit Zotero or copy DB if locked)."
