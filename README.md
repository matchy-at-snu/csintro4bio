# Introduction to Computer Science for Biologists / 생물학자를 위한 전산학 개론

Course taken at Seoul National University in 2019 spring semester.

## Disclaimers

1. The author did not intentionally make the choice of using [Hungarian notation](https://en.wikipedia.org/wiki/Hungarian_notation) + small camel case to name variables. The author is very well aware that snake case is the canonical way to name variables in Python. It was a requirement forced by the lecturer.
2. Code reusing was not allowed by the lecturer as the submission was required to be a single file. The author is aware that some classes were defined multiple times.
3. The author got B+ in this course. It doesn't mean the code produces wrong results, though. That said, don't rely on this repo to get a good grade.

## How to run

0. Clone the repo
1. Set up proper environemt either via:
   1.1. `pip install -r requirements.txt`
   1.2. `conda env create -f environment.yml`
   1.3. Manually install the required packages, mission 4 requires `scipy`, mission 5 requires `statsmodels`. That's it.
2. Run `download_data.sh` to download the data
3. `python mission/00/Mission0.py` to run the code. Replace `00` and `Mission0` with the mission number you want to run.
   2.1. Alternatively, use `run_all.sh` to run all the missions at once
4. See stdout or check "mission/xx/solutions" for the output.
