@echo off
cd %cd%
docker run -it --rm -v $(pwd):/data p_tree:latest /data/pipeline.sh
pause