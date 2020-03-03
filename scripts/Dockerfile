FROM python:3.8

WORKDIR /usr/src/app

COPY requirements.txt ./

# RUN apk add --no-cache build-base lapack-dev

RUN pip search pytables
RUN pip install --no-cache-dir -r requirements.txt

COPY analysis.py .
COPY analysis_helper_functions.py .


ENTRYPOINT ["python", "./analysis.py"]