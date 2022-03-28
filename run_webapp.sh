#!/bin/bash

streamlit run webapp.py \
    --server.headless true \
    --browser.gatherUsageStats false \
    --browser.serverAddress "127.0.0.1"

