
FROM python:3.9
# FROM ubuntu:latest

# RUN apt-get update && apt-get install -y \
#     python3.10 \ 
#     python3-pip 

WORKDIR /code

# # RUN apk --no-cache add musl-dev linux-headers g++
# COPY requirements.txt requirements.txt

# COPY . /code 

# FROM python:3.9
# WORKDIR /code 

# # RUN apk --no-cache add musl-dev linux-headers g++
# COPY requirements.txt requirements.txt

# RUN ls
# CMD [ "python3", "app.py" ] 

# Remove SSH keys
RUN rm -rf /root/.ssh/

# Add the rest of the files
ADD . .

# CMD python manage.py runserver