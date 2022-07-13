#!/bin/bash

# Launch the streamlit-based webapp!

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd ${SCRIPT_DIR}

streamlit run Welcome.py \
    --server.headless true \
    --browser.gatherUsageStats false \
    --browser.serverAddress "127.0.0.1" \
    Welcome.py
