@echo off
mkdir input && mkdir output && mkdir reference
docker build -t p_tree:latest .
echo Docker build successful!
pause