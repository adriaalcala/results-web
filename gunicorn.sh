#!/bin/bash

[[ "$APP_ENV" = "production" ]] && WORKERS="4" || WORKERS="1"
[[ "$APP_ENV" = "production" ]] && RELOAD="" || RELOAD="--reload"

/usr/local/bin/gunicorn app.server:app \
    --chdir=/opt \
    -w ${WORKERS} \
    -b 0.0.0.0:8080 \
    --worker-class aiohttp.GunicornUVLoopWebWorker \
    ${RELOAD} \
    --access-logfile "-" \
    --access-logformat '%a %t "%r" %s %b "%{Referer}i" "%{User-Agent}i"'