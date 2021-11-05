# Use ubuntu as a base image
FROM ubuntu:20.04

# Label container as authored by Jacob Oakman
LABEL org.opencontainers.image.authors="Jacob Oakman"

# Update and Upgrade apt packages
RUN apt update -y
RUN apt upgrade -y

# Install additional packages for development
RUN DEBIAN_FRONTEND=noninteractive apt install -y git man curl build-essential gdb python3 libspdlog-dev

# Change the working directory
WORKDIR "/root"
