@echo off
cd %cd%
docker run -it --rm -v D:\Git_Project\Tree:/data p_tree:latest /data/pipeline.sh
pause