#!/bin/bash

# Launch the streamlit-based webapp, either locally on your computer or automatically
# on binder via https://gitlab.mpcdf.mpg.de/khr/jupyter-streamlit-proxy

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd ${SCRIPT_DIR}

streamlit run Welcome.py \
    --server.headless true \
    --browser.gatherUsageStats false \
    --browser.serverAddress "0.0.0.0" \
    "$@"

#    --browser.serverAddress "127.0.0.1" "$@"

