
FROM python:3.9
WORKDIR /code
ENV FLASK_APP app.py
ENV FLASK_RUN_HOST 0.0.0.0
# RUN apk add --no-cache gcc musl-dev linux-headers
COPY requirements.txt requirements.txt

RUN pip install --upgrade pip
RUN pip install --no-cache-dir -r requirements.txt 

COPY . /code
CMD ["flask", "run"]


# FROM ubuntu:latest

# RUN apt-get update && apt-get install -y \
#     python3.10 \ 
#     python3-pip 

# FROM python:3.9
# WORKDIR /code 

# # RUN apk --no-cache add musl-dev linux-headers g++
# COPY requirements.txt requirements.txt

# COPY . /code 

# RUN ls
# CMD [ "python3", "app.py" ] 
