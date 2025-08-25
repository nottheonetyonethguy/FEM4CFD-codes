# official texlive image
FROM texlive/texlive:latest

# dependencies
RUN apt-get update && \
        apt-get install -y --no-install-recommends \
        zathura \
        zathura-pdf-poppler \
        neovim \
        git \
        curl \
        wget \
        xdg-utils \
        dbus-x11 \
        nodejs \
        npm \
	build-essential \
	libwayland-client0 \
	libwayland-cursor0 \
	libxkbcommon0 \
        && rm -rf /var/lib/apt/lists*

# create a user
RUN useradd -m nottheonetyonethguy && \
        mkdir -p /home/nottheonetyonethguy/.config/nvim && \
        chown -R nottheonetyonethguy:nottheonetyonethguy /home/nottheonetyonethguy

RUN npm install -g tree-sitter-cli

# switch user
USER nottheonetyonethguy
WORKDIR /home/nottheonetyonethguy

# wayland specific 
ENV XDG_RUNTIME_DIR=/tmp/runtime-nottheonetyonethguy
ENV WAYLAND_DISPLAY=wayland-1
ENV QT_QPA_PLATFORM=wayland
ENV GDK_BACKEND=wayland
ENV CLUTTER_BACKEND=wayland
ENV SDL_VIDEODRIVER=wayland
ENV MOZ_ENABLE_WAYLAND=1

# create runtime
RUN mkdir -p ${XDG_RUNTIME_DIR} && \
    chmod 700 ${XDG_RUNTIME_DIR}

# env
ENV DISPLAY=host.docker.internal:0
ENV DBUS_SESSION_BUS_ADDRESS=unix:path=/run/user/1000/bus

CMD [ "bash" ]
