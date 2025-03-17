@echo off
cd D:\Git_Project\Tree
mkdir input output 2>nul
docker run -it --rm -v D:\Git_Project\Tree:/data p_tree:latest /data/pipeline.sh
pause