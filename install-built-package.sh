#/bin/bash
# https://packaging.python.org/distributing/#working-in-development-mode

#only works after running build-package.sh

#python3 -mpip install opencv-python
python3 -mpip install pandas jsonlib-python3 datetime

pip3 install --no-cache-dir ./package-project/src/dist/pytec-1.0-py2.py3-none-any.whl
