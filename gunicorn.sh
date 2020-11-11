#!/bin/bash

/usr/local/bin/gunicorn app:server \
    --chdir=/opt \
    -b 0.0.0.0:8080 \
    --reload
