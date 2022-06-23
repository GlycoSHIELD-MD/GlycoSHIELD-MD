#!/bin/bash

streamlit run Welcome.py \
    --server.headless true \
    --browser.gatherUsageStats false \
    --browser.serverAddress "127.0.0.1"

