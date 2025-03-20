@echo off
docker run -it --rm -v %cd%:/data p_tree:latest /data/pipeline.sh
pause