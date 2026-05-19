#!/usr/bin/env bash
set -euo pipefail

APP_PORT="${APP_PORT:-8000}"
DATASETTE_PORT="${DATASETTE_PORT:-8001}"
DB="${DB:-pipeline/ptm_disease.db}"
HOST="${HOST:-127.0.0.1}"

cd "$(dirname "$0")"

if [[ ! -f "$DB" ]]; then
  echo "Database not found at $DB" >&2
  exit 1
fi

cleanup() {
  echo
  echo "[serve] stopping..."
  kill 0 2>/dev/null || true
}
trap cleanup EXIT INT TERM

echo "[serve] starting FastAPI on http://$HOST:$APP_PORT (dashboard + /api/*)"
uv run uvicorn server.app:app --host "$HOST" --port "$APP_PORT" 2>&1 | sed 's/^/[app] /' &

echo "[serve] starting datasette on http://$HOST:$DATASETTE_PORT (raw SQL access)"
uv run datasette "$DB" --host "$HOST" --port "$DATASETTE_PORT" 2>&1 | sed 's/^/[datasette] /' &

sleep 1

echo "[serve] starting cloudflared tunnel -> http://$HOST:$APP_PORT"
cloudflared tunnel --url "http://$HOST:$APP_PORT" 2>&1 | sed 's/^/[cloudflared] /' &

wait
