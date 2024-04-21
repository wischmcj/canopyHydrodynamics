
FROM python:3.9
# FROM ubuntu:latest

# RUN apt-get update && apt-get install -y \
#     python3.10 \ 
#     python3-pip 


# FROM python:3.6-slim

ARG ssh_prv_key
ARG ssh_pub_key

RUN apt-get update && \
    apt-get install -y \
        git \
        openssh-server \
        libmysqlclient-dev

# Authorize SSH Host
RUN mkdir -p /root/.ssh && \
    chmod 0700 /root/.ssh
# See: https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/githubs-ssh-key-fingerprints
COPY known_hosts > /root/.ssh/known_hosts

# Add the keys and set permissions
RUN echo "$ssh_prv_key" > /root/.ssh/id_rsa && \
    echo "$ssh_pub_key" > /root/.ssh/id_rsa.pub && \
    chmod 600 /root/.ssh/id_rsa && \
    chmod 600 /root/.ssh/id_rsa.pub

# Avoid cache purge by adding requirements first
# ADD ./requirements.txt /app/requirements.txt
COPY requirements.txt requirements.txt

WORKDIR /code


# RUN  pip install --upgrade pip
RUN pip install -r requirements.txt

# Remove SSH keys
RUN rm -rf /root/.ssh/

# Add the rest of the files
ADD . .

COPY . /code   
CMD [ "python3", "app.py" ] 
# CMD python manage.py runserver