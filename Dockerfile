
FROM python:3.9
# FROM ubuntu:latest

# RUN apt-get update && apt-get install -y \
#     python3.10 \ 
#     python3-pip 

WORKDIR /code


# RUN apk --no-cache add musl-dev linux-headers g++
COPY requirements.txt requirements.txt

# RUN  pip install --upgrade pip
RUN pip install -r requirements.txt

COPY . /code   
CMD [ "python3", "app.py" ] 