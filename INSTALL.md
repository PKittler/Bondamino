# Installing and running Bondamino

This guide explains how to install Bondamino, how to configure it and how to keep it running:
updates, backups, security and troubleshooting. It is written for people who have never
installed server software before, but it is also meant to be complete enough for system
administrators. If you only want to know what Bondamino does, start with the
[README](README.md); how to use the application is described in the [tutorial](TUTORIAL.md)
and in the [documentation](DOCUMENTATION.md).

Bondamino is a web application written in Python with the Django framework. You install it once
on a computer (your own laptop or a server), and then everybody uses it in the web browser.

## Contents

- [Before you start](#before-you-start)
  - [Basic terms](#basic-terms)
  - [Opening a terminal](#opening-a-terminal)
  - [Installing the tools](#installing-the-tools)
- [Choosing an installation method](#choosing-an-installation-method)
- [Choosing a database](#choosing-a-database)
- [Quick try-out on your own computer](#quick-try-out-on-your-own-computer)
  - [Try-out with Python](#try-out-with-python)
  - [Try-out with Docker](#try-out-with-docker)
- [Standalone installation on a Linux server](#standalone-installation-on-a-linux-server)
- [The guided setup in detail](#the-guided-setup-in-detail)
- [Standalone installation on Windows](#standalone-installation-on-windows)
- [Reverse proxy and HTTPS](#reverse-proxy-and-https)
- [Docker Compose](#docker-compose)
- [Portainer](#portainer)
- [Other container platforms](#other-container-platforms)
- [After the installation](#after-the-installation)
- [Configuration reference](#configuration-reference)
- [Updating](#updating)
- [Backup and restore](#backup-and-restore)
- [Moving to another server](#moving-to-another-server)
- [Uninstalling](#uninstalling)
- [Security notes for operators](#security-notes-for-operators)
- [Troubleshooting](#troubleshooting)

---

## Before you start

### Basic terms

You will meet these terms in this guide. Each is explained briefly here; you do not need to
know more than this to follow the steps.

| Term | Meaning |
|---|---|
| Terminal | A window in which you type commands instead of clicking. Also called command line, shell, PowerShell or console. |
| Server | A computer that is always on and that other people reach over the network. It can be a physical machine or a virtual machine (VM) provided by your IT department. |
| Python | The programming language Bondamino is written in. You need it for the [standalone installation](#standalone-installation-on-a-linux-server), not for Docker. |
| Virtual environment | A folder (`.venv`) that holds the Python packages of one application, separate from the rest of the computer. |
| Git | A program that downloads the source code and later fetches updates. |
| Database server | A program (PostgreSQL or MariaDB) that stores the data (users, groups, files, notes) and answers queries from Bondamino. SQLite is a database in a single file and needs no separate server. |
| Web server (WSGI server) | The program that runs Bondamino and answers browser requests: gunicorn on Linux/macOS, waitress on Windows. |
| Reverse proxy | A web server in front of Bondamino (nginx, Apache, Caddy, IIS, Traefik). It receives the requests from the internet, handles HTTPS and passes them on to Bondamino. |
| HTTPS, certificate | HTTPS encrypts the connection between browser and server. It needs a certificate for your domain, e.g. free from Let's Encrypt or from your institute. |
| Domain, host name | The address under which Bondamino is reached, e.g. `bondamino.example.org`. |
| Port | A number that selects a program on a computer. Bondamino listens on port 8000; browsers use port 443 for HTTPS and 80 for HTTP. |
| Environment variable | A named setting (e.g. `ALLOWED_HOSTS=bondamino.example.org`) that is given to a program from outside. Bondamino reads its whole configuration from environment variables or from a file with such lines (`.env`). |
| Docker, container, image | Docker runs programs in isolated "containers". An image is the packaged program from which a container is started. Docker Compose starts several containers (e.g. Bondamino and its database) from one file. |
| Volume | A storage area that Docker keeps when a container is replaced, e.g. during an update. Bondamino keeps its uploaded files there. |
| Portainer | A web interface for managing Docker. |

### Opening a terminal

- **Windows:** open the Start menu, type `PowerShell` and click **Windows PowerShell** (or
  **Terminal**). For commands marked "as administrator", right-click it and choose **Run as
  administrator**.
- **macOS:** press `Cmd` + `Space`, type `Terminal` and press `Enter`.
- **Linux desktop:** press `Ctrl` + `Alt` + `T`, or find **Terminal** in the application menu.
- **Linux server:** connect from your own computer with SSH (works the same in PowerShell, the
  macOS Terminal and Linux):

  ```sh
  ssh your-user-name@server.example.org
  ```

Type each command exactly as shown and press `Enter`. Lines in the code blocks that start with
`#` are comments; you do not need to type them. Commands that start with `sudo` run with
administrator rights on Linux and ask for your password; while you type a password, nothing is
shown on the screen, not even dots. This is normal.

### Installing the tools

Which tools you need depends on the installation method:

| Method | Python 3.12 or newer | Git | Docker |
|---|---|---|---|
| Try-out with Python, standalone installation | yes | recommended | no |
| Try-out with Docker, Docker Compose, Portainer | no | recommended | yes |

#### Python

Bondamino needs **Python 3.12 or newer**.

- **Windows:** download the installer from <https://www.python.org/downloads/> and, on the
  first screen, tick **Add python.exe to PATH**. Alternatively, in PowerShell:
  `winget install Python.Python.3.13`.
- **macOS:** download the installer from <https://www.python.org/downloads/> (or, if you use
  Homebrew, `brew install python@3.13`).
- **Linux (Debian/Ubuntu):**

  ```sh
  sudo apt update
  sudo apt install python3 python3-venv python3-pip
  ```

  Ubuntu 24.04 and Debian 13 ship a suitable version. Debian 12 ships Python 3.11, which is too
  old; use a newer system or [Docker](#docker-compose).

Check the version (close and reopen the terminal after installing):

```sh
python3 --version     # macOS and Linux
py --version          # Windows
```

Expected output: `Python 3.12.x` or higher.

#### Git

- **Windows:** download from <https://git-scm.com/download/win> or run
  `winget install --id Git.Git -e`.
- **macOS:** type `git --version` in the Terminal; if Git is missing, macOS offers to install
  the "command line developer tools". Accept.
- **Linux:** `sudo apt install git`

Check with `git --version` (expected: `git version 2.x`).

**Without Git:** on <https://gitlab.com/d11m-group/bondamino> click **Code**, then under
**Download source code** choose **zip**, and extract the archive. The extracted folder is called
`bondamino-main`; rename it to `Bondamino` or use its name wherever this guide says `Bondamino`. Updates are easier with Git, because with a ZIP you have to download and copy
the new version yourself (see [Updating](#updating)).

#### Docker

- **Windows and macOS:** install **Docker Desktop** from
  <https://www.docker.com/products/docker-desktop/> and start it once (on Windows it uses WSL 2
  and may ask you to restart the computer). Docker Desktop is free for personal use, education
  and small organisations; larger organisations may need a paid subscription, so check the
  current terms with your IT department.
- **Linux server:** install **Docker Engine** following
  <https://docs.docker.com/engine/install/> (choose your distribution). Then allow your user to
  run Docker without `sudo` and log out and in again:

  ```sh
  sudo usermod -aG docker $USER
  ```

Check that Docker works:

```sh
docker --version
docker compose version
docker run --rm hello-world
```

Expected: `Docker version 24` or newer, `Docker Compose version v2.x`, and a message
"Hello from Docker!". If `docker compose` is unknown, you have an old Compose v1 installation
(`docker-compose`); install the Compose v2 plugin.

---

## Choosing an installation method

| Your situation | Recommended method |
|---|---|
| I want to try Bondamino on my own computer | [Try-out with Docker](#try-out-with-docker) or [with Python](#try-out-with-python) and SQLite |
| A small working group, one server, nobody knows Docker | [Standalone on Linux](#standalone-installation-on-a-linux-server) with PostgreSQL (or SQLite), or [on Windows](#standalone-installation-on-windows) |
| An institute server run by an IT department | [Docker Compose](#docker-compose) with PostgreSQL, or standalone with PostgreSQL, behind the existing reverse proxy |
| There is already a Docker host or Portainer | [Docker Compose](#docker-compose) or [Portainer](#portainer) |
| Kubernetes, Nomad, Synology, Unraid, cloud services | [Other container platforms](#other-container-platforms) |

All methods are configured the same way, through environment variables (see the
[configuration reference](#configuration-reference)). You can switch methods later by moving the
data (see [Moving to another server](#moving-to-another-server)).

**About the prebuilt image:** the repository contains a build pipeline (`.gitlab-ci.yml`) that
publishes a ready container image as `registry.gitlab.com/d11m-group/bondamino`. Whether you can
download it depends on whether the maintainer has made the project's container registry public. This guide therefore builds the image from the source code on your machine, which
works in every case. If the maintainer has published the image, you can use it instead; this is
mentioned at the relevant places.

**Requirements for a server:** a small virtual machine (2 CPU cores, 2 GB RAM) is enough for a
working group. Disk space depends on the uploaded PDB files (up to 20 MB each). The server needs
internet access during installation and updates (packages, images), and, if structures are
loaded by PDB ID, access to `https://files.rcsb.org`. Sending invitations needs a mail server
(SMTP); ask your IT department for its address and credentials.

---

## Choosing a database

| Database | Use |
|---|---|
| **PostgreSQL** (version 14 or newer) | Recommended for production. Robust with many simultaneous users; the Docker Compose default. |
| **MariaDB** (10.5 or newer) / **MySQL** (8.0 or newer) | Fully supported. Choose it if your institute already runs MariaDB/MySQL. |
| **SQLite** | Testing, single users and small groups. A single file, no database server needed. Writes are handled one after another, so it is not meant for many concurrent users. |

MongoDB and DuckDB are not supported: Bondamino's data (users, groups, files, notes) is relational
and depends on the guarantees of a relational database; MongoDB is a document store that Django
does not support, and DuckDB is designed for analytics, not for the many small concurrent writes
of a web application.

The database driver for PostgreSQL (`psycopg`) is licensed under the LGPL and the driver for
MariaDB/MySQL (`mysqlclient`) under the GPL; the Docker image contains both. See
[LICENSE.md](LICENSE.md) for what this means.

---

## Quick try-out on your own computer

> **Not for production.** The try-out runs Bondamino on your own computer for testing, without
> HTTPS and with emails only printed to the terminal. Do not make it reachable from the
> internet and do not store real data in it.

### Try-out with Python

You need [Python](#python) and [Git](#git).

1. Open a [terminal](#opening-a-terminal) and download Bondamino:

   ```sh
   git clone https://gitlab.com/d11m-group/bondamino.git Bondamino
   cd Bondamino
   ```

   Expected: `Cloning into 'Bondamino'...` followed by `done`. The folder `Bondamino` is
   created in the current folder (usually your home folder).

2. Create and activate a virtual environment.

   macOS and Linux:

   ```sh
   python3 -m venv .venv
   . .venv/bin/activate
   ```

   Windows (PowerShell):

   ```powershell
   py -m venv .venv
   .venv\Scripts\Activate.ps1
   ```

   The prompt now starts with `(.venv)`. Repeat the activation each time you open a new terminal.
   If PowerShell reports that "running scripts is disabled on this system", run
   `Set-ExecutionPolicy -Scope CurrentUser RemoteSigned`, confirm with `Y` and try again.

3. Install the packages (SQLite needs no extra driver):

   ```sh
   pip install -r requirements.txt
   ```

   Expected at the end: `Successfully installed Django-5.2... tzdata-... whitenoise-...`.
   (`tzdata` is the time zone database; Windows has none of its own.)

4. Run the [guided setup](#the-guided-setup-in-detail):

   ```sh
   cd bondamino
   python manage.py setup
   ```

   Sensible answers for a try-out (press `Enter` to accept the value in brackets):

   | Question | Answer |
   |---|---|
   | Database | `1` (SQLite) |
   | Database file | `Enter` (`db.sqlite3`) |
   | Domain(s) | `Enter` (`localhost`; `127.0.0.1` is added automatically) |
   | HTTPS | `n` |
   | Development mode (DEBUG) | `y` |
   | Default time zone | `Enter` (`UTC`), or e.g. `Europe/London` |
   | Sending (email) | `2` (print emails in the console) |
   | Write ...? | `y` |

   The setup then creates the database and asks for the first administrator: a user name, an
   email address and a password (typed twice, invisible). Expected at the end:
   `Bondamino is set up.`

5. Start the development server:

   ```sh
   python manage.py runserver
   ```

   Expected: `Starting development server at http://127.0.0.1:8000/`. Open
   <http://localhost:8000/> in your browser and sign in with the administrator. Stop the server
   with `Ctrl` + `C`.

Next time: open a terminal, `cd Bondamino`, activate the virtual environment (step 2, second
line), `cd bondamino`, `python manage.py runserver`.

### Try-out with Docker

You need [Docker](#docker) and, recommended, [Git](#git).

1. Download Bondamino:

   ```sh
   git clone https://gitlab.com/d11m-group/bondamino.git Bondamino
   cd Bondamino
   ```

2. Create the configuration file `.env` from the template.

   macOS and Linux: `cp .env.example .env`, then edit it, e.g. with `nano .env`
   (save with `Ctrl` + `O`, `Enter`, quit with `Ctrl` + `X`) or `open -e .env` on macOS.

   Windows: `Copy-Item .env.example .env`, then `notepad .env`.

   Change these lines (leave everything else as it is):

   ```ini
   ADMIN_USERNAME=admin
   ADMIN_EMAIL=you@example.org
   ADMIN_PASSWORD=a-long-password-of-your-choice
   ```

   `EMAIL_BACKEND` stays empty, so Bondamino prints emails into its log instead of sending them.

3. Build the image and start Bondamino with SQLite:

   ```sh
   docker compose -f docker-compose.sqlite.yml up -d --build
   ```

   The first build downloads and compiles the dependencies and takes a few minutes. You may see
   a message that the image `registry.gitlab.com/d11m-group/bondamino` could not be pulled; this is expected,
   the image is built locally instead. Expected at the end: `Container bondamino-bondamino-1
   Started`.

4. Watch the start-up (stop watching with `Ctrl` + `C`; Bondamino keeps running):

   ```sh
   docker compose -f docker-compose.sqlite.yml logs -f bondamino
   ```

   Expected: `Generated a new secret key in /data/secret_key`, a list of `Applying ...` lines,
   `Superuser created successfully.` and `Listening at: http://0.0.0.0:8000`.

5. Open <http://localhost:8000/> and sign in with `ADMIN_USERNAME` and `ADMIN_PASSWORD`.

Stop it with `docker compose -f docker-compose.sqlite.yml down` (the data is kept) and start it
again with `docker compose -f docker-compose.sqlite.yml up -d`. To remove it completely, see
[Uninstalling](#uninstalling).

---

## Standalone installation on a Linux server

This section installs Bondamino directly on a Debian or Ubuntu server: gunicorn as web server,
systemd to start it automatically, nginx as reverse proxy and a Let's Encrypt certificate for
HTTPS. Other distributions work the same way with their own package manager (`dnf` instead of
`apt`, etc.).

The examples use these names; replace them with your own:

| Example | Meaning |
|---|---|
| `bondamino.example.org` | The domain under which Bondamino will be reached. It must point to the server (ask your IT department for a DNS entry). |
| `/opt/Bondamino` | Installation folder |
| `bondamino` | System user that runs Bondamino, and name of the database and database user |

### Step 1: Prepare the server

Connect to the server ([SSH](#opening-a-terminal)) and install the required packages:

```sh
sudo apt update
sudo apt upgrade
sudo apt install python3 python3-venv python3-pip git nginx
python3 --version
```

The last command must show Python 3.12 or newer.

Create a system user that runs Bondamino (it cannot log in, which limits the damage if something
goes wrong):

```sh
sudo useradd --system --create-home --home-dir /home/bondamino --shell /usr/sbin/nologin bondamino
```

### Step 2: Create the database

Skip this step for SQLite.

#### PostgreSQL

```sh
sudo apt install postgresql
sudo -u postgres psql
```

The prompt changes to `postgres=#`. Choose a password (letters and digits avoid quoting problems)
and type:

```sql
CREATE USER bondamino WITH PASSWORD 'choose-a-password';
CREATE DATABASE bondamino OWNER bondamino;
\q
```

Expected: `CREATE ROLE` and `CREATE DATABASE`. Check that the new user can sign in (it asks for
the password; expected output is a small table with `1`):

```sh
psql -h localhost -U bondamino -d bondamino -c 'SELECT 1;'
```

#### MariaDB

```sh
sudo apt install mariadb-server
sudo mariadb
```

The prompt changes to `MariaDB [(none)]>`. Type:

```sql
CREATE DATABASE bondamino CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;
CREATE USER 'bondamino'@'localhost' IDENTIFIED BY 'choose-a-password';
GRANT ALL PRIVILEGES ON bondamino.* TO 'bondamino'@'localhost';
EXIT;
```

Expected: `Query OK` after each statement. Check:

```sh
mariadb -u bondamino -p bondamino -e 'SELECT 1;'
```

The MariaDB/MySQL driver is compiled during installation and needs a compiler and the client
headers:

```sh
sudo apt install build-essential pkg-config python3-dev default-libmysqlclient-dev
```

#### A database on another server

If the database runs on a separate server (e.g. a central database service of your institute),
ask for host name, port, database name, user and password, and whether TLS is required. You will
enter these in the setup. For PostgreSQL with TLS, add `?sslmode=require` to `DATABASE_URL`
afterwards (see [DATABASE_URL](#database_url)).

### Step 3: Download Bondamino and install the packages

```sh
sudo git clone https://gitlab.com/d11m-group/bondamino.git /opt/Bondamino
sudo chown -R bondamino:bondamino /opt/Bondamino
```

From now on, work as the `bondamino` user so that all files belong to it:

```sh
sudo -u bondamino -H bash
cd /opt/Bondamino
python3 -m venv .venv
. .venv/bin/activate
```

The prompt now starts with `(.venv)`. Install the packages for your database (only one of the
three lines):

```sh
pip install -r requirements-postgresql.txt    # PostgreSQL
pip install -r requirements-mariadb.txt       # MariaDB / MySQL
pip install -r requirements.txt               # SQLite
```

Expected at the end: `Successfully installed ...` including `Django`, `gunicorn` and the driver
(`psycopg` or `mysqlclient`). If `mysqlclient` fails, see
[mysqlclient cannot be installed](#mysqlclient-cannot-be-installed).

### Step 4: Run the guided setup

Still as the `bondamino` user with the virtual environment active:

```sh
cd /opt/Bondamino/bondamino
python manage.py setup
```

For a server, answer:

| Question | Answer |
|---|---|
| Database | `2` (PostgreSQL) or `3` (MariaDB), then server `localhost`, the port, name, user and password from step 2 |
| Domain(s) | `bondamino.example.org` |
| HTTPS | `y` (nginx will provide HTTPS in step 7) |
| Development mode (DEBUG) | `n` |
| Default time zone | e.g. `Europe/Berlin` (or `Enter` for `UTC`) |
| Email | `1` (SMTP) with the data of your mail server, and send a test email to yourself |

All questions are explained in [The guided setup in detail](#the-guided-setup-in-detail). The
result is the file `/opt/Bondamino/bondamino/.env`, the database tables, the collected static
files (styles, scripts) and the first administrator.

### Step 5: Test the web server

Start gunicorn by hand to check that everything works:

```sh
gunicorn bondamino.wsgi:application --bind 127.0.0.1:8000 --workers 3 --timeout 120
```

Expected: `Listening at: http://127.0.0.1:8000` and `Booting worker with pid: ...` (three
times). In a second SSH session, run `curl -I -H 'Host: bondamino.example.org' http://127.0.0.1:8000/login/`;
expected first line: `HTTP/1.1 200 OK`. Stop gunicorn with `Ctrl` + `C`, then leave the
`bondamino` user with `exit`.

`--bind 127.0.0.1:8000` makes Bondamino reachable only from the server itself; the outside world
reaches it through nginx. `--workers` is the number of processes (about two per CPU core is
common); `--timeout 120` allows slow requests (e.g. large files) up to 120 seconds.

### Step 6: Start Bondamino automatically with systemd

systemd is the Linux service manager; it starts Bondamino at boot and restarts it after a crash.
Create the service file:

```sh
sudo nano /etc/systemd/system/bondamino.service
```

with this content:

```ini
[Unit]
Description=Bondamino
After=network.target postgresql.service mariadb.service

[Service]
User=bondamino
Group=bondamino
WorkingDirectory=/opt/Bondamino/bondamino
ExecStart=/opt/Bondamino/.venv/bin/gunicorn bondamino.wsgi:application --bind 127.0.0.1:8000 --workers 3 --timeout 120 --access-logfile -
Restart=on-failure

[Install]
WantedBy=multi-user.target
```

Enable and start it, then check its status:

```sh
sudo systemctl daemon-reload
sudo systemctl enable --now bondamino
sudo systemctl status bondamino
```

Expected: `Active: active (running)`. Leave the status view with `q`. The log (including emails,
if they are only printed to the console) is shown with:

```sh
sudo journalctl -u bondamino -f
```

### Step 7: nginx and HTTPS

nginx receives the requests on ports 80 and 443 and passes them on to gunicorn. Create the site
configuration:

```sh
sudo nano /etc/nginx/sites-available/bondamino
```

```nginx
server {
    listen 80;
    server_name bondamino.example.org;

    # Uploads are limited to 20 MB by Bondamino; allow a little more for the form data
    client_max_body_size 25m;

    location / {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_read_timeout 120s;
    }
}
```

Enable it and check the configuration:

```sh
sudo ln -s /etc/nginx/sites-available/bondamino /etc/nginx/sites-enabled/bondamino
sudo nginx -t
sudo systemctl reload nginx
```

Expected from `nginx -t`: `syntax is ok` and `test is successful`.

Now get a free certificate from Let's Encrypt with certbot. This only works if the server is
reachable from the internet on port 80 under `bondamino.example.org`:

```sh
sudo apt install certbot python3-certbot-nginx
sudo certbot --nginx -d bondamino.example.org
```

certbot asks for an email address (for expiry warnings) and adds the HTTPS part (`listen 443
ssl`, certificate paths) and a redirect from HTTP to HTTPS to the file above. It also renews the
certificate automatically; test that with `sudo certbot renew --dry-run`.

**Server only reachable inside the institute network:** Let's Encrypt cannot reach it. Ask your
IT department for a certificate and add it yourself:

```nginx
server {
    listen 443 ssl;
    server_name bondamino.example.org;
    ssl_certificate     /etc/ssl/certs/bondamino.crt;      # certificate incl. intermediate certificates
    ssl_certificate_key /etc/ssl/private/bondamino.key;
    client_max_body_size 25m;
    location / {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Forwarded-Proto $scheme;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_read_timeout 120s;
    }
}
server {
    listen 80;
    server_name bondamino.example.org;
    return 301 https://$host$request_uri;
}
```

Recommended: enable HSTS (browsers then always use HTTPS for this domain) by adding this line
inside the `listen 443` server block, then `sudo nginx -t && sudo systemctl reload nginx`:

```nginx
    add_header Strict-Transport-Security "max-age=31536000" always;
```

### Step 8: Firewall

If the server uses the `ufw` firewall, allow SSH, HTTP and HTTPS. Port 8000 stays closed; it is
only used inside the server.

```sh
sudo ufw allow OpenSSH
sudo ufw allow 'Nginx Full'
sudo ufw enable
```

Open `https://bondamino.example.org/` in your browser. Continue with
[After the installation](#after-the-installation).

Apache and Caddy instead of nginx are described in [Reverse proxy and HTTPS](#reverse-proxy-and-https).

---

## The guided setup in detail

`python manage.py setup` (run in the folder `bondamino`, the one that contains `manage.py`)
asks for the configuration step by step, tests it, and writes the file `bondamino/.env`. It is
only needed for the standalone installation; Docker and Portainer are configured with
environment variables instead.

Press `Enter` to accept the value shown in brackets. Yes/no questions accept `y`/`n` (and
`j`/`ja`/`nein`). Passwords are not shown while you type.

**1. Database**

| Question | Meaning |
|---|---|
| Database `1`/`2`/`3` | SQLite, PostgreSQL or MariaDB/MySQL (see [Choosing a database](#choosing-a-database)). |
| Database file (SQLite) | File name relative to the folder of `manage.py` (`db.sqlite3`) or an absolute path. The folder must exist. |
| Server | Host name of the database server; `localhost` if it runs on the same machine. |
| Port | `5432` (PostgreSQL) or `3306` (MariaDB) unless your database uses another one. |
| Database name, Database user, Database password | As created in [step 2](#step-2-create-the-database). |

The setup then connects to the database. `Connection OK.` means everything is right. On
`Connection failed: ...` it shows the error and offers to enter the settings again. If the driver
is missing, it tells you which requirements file to install (from the repository root). Special
characters in user and password are encoded automatically.

**2. Web address**

| Question | Meaning |
|---|---|
| Domain(s) under which Bondamino is reached | The host names that appear in the browser's address bar, comma-separated, e.g. `bondamino.example.org`. Requests for any other name are rejected with "400 Bad Request". For a try-out: `localhost`. If the list contains `localhost`, the setup adds `127.0.0.1` automatically. |
| Is Bondamino served via HTTPS (through a reverse proxy)? | `y` if a reverse proxy provides HTTPS (the normal case on a server). Cookies are then only sent over HTTPS, and the proxy's information about protocol and client address is trusted. `n` for a try-out on `http://localhost`. |
| Development mode (DEBUG)? | `n` on every server. `y` shows detailed error pages with internal information; only for a try-out or development. |
| Default time zone | Used to display dates and times, e.g. `UTC` (the default), `Europe/Berlin`, `Europe/London`, `America/New_York` (names from the IANA time zone database). An unknown name is rejected and the question is asked again. Each user can choose their own time zone in their account settings; this is the default for everybody who has not. |

**3. Email (invitations, forgotten passwords)**

| Question | Meaning |
|---|---|
| Sending `1`/`2` | `1` sends emails via an SMTP server. `2` only prints them to the server console (try-out). Without real emails, users cannot receive invitations or reset their password. |
| SMTP server | Host name of the mail server, e.g. `smtp.example.org` (from your IT department or mail provider). |
| Encryption | `1` STARTTLS (usually port 587), `2` SSL/TLS (usually port 465), `3` none (usually port 25, only inside trusted networks). |
| Port | Suggested to match the encryption. |
| SMTP user, SMTP password | Login of the mailbox used for sending; leave the user empty if the server accepts mail without login. |
| Sender address | The "From" of all emails, e.g. `Bondamino <bondamino@example.org>`. Use an address your mail server is allowed to send from; many servers reject other senders. |
| Send a test email to | Your own address. `Test email sent to ...` means the settings work; otherwise the error is shown. |

**4. Writing the file**

The setup creates a random `SECRET_KEY` (only the first time), shows a summary with the passwords
masked and asks `Write .../.env?`. The previous file is kept as `.env.bak`. Values that you added
to `.env` by hand are kept. On Linux and macOS the file is made readable for its owner only.

**5. Preparing the database**

Finally the setup creates or updates the database tables, creates the cache table used for
limiting sign-in attempts, collects the static files and, if there is no administrator yet, asks
for the first one (user name, email address, password twice; the prompts may partly appear in
German). A password that is too short or too common triggers a warning; answer `N` and choose a
better one. Expected at the end: `Bondamino is set up.`

**Running it again:** you can run the setup at any time to change the configuration; your
current values are offered as defaults. Restart the web server afterwards (e.g.
`sudo systemctl restart bondamino`). Options:

| Command | Effect |
|---|---|
| `python manage.py setup` | Full guided setup |
| `python manage.py setup --skip-finish` | Only write `.env`, do not touch the database |
| `python manage.py setup --finish` | Only prepare the database with the current configuration (migrations, cache table, static files, first administrator) |
| `python manage.py setup --env-file /path/to/file` | Write another file instead of `bondamino/.env` (Bondamino reads it if `BONDAMINO_ENV_FILE` points to it) |

You can also edit `.env` with a text editor; the meaning of every line is in the
[configuration reference](#configuration-reference).

---

## Standalone installation on Windows

Bondamino runs on Windows (e.g. Windows Server or a Windows PC in the lab) with the waitress web
server, because gunicorn does not run on Windows. The examples install it to `C:\Bondamino`.
Use PowerShell **as administrator** for the commands in this section.

1. Install [Python](#python) (tick **Add python.exe to PATH**) and [Git](#git).

2. Install a database server if you do not use SQLite:
   - PostgreSQL: installer from <https://www.postgresql.org/download/windows/>. Remember the
     password of the `postgres` user. Then open **SQL Shell (psql)** from the Start menu,
     accept the defaults with `Enter`, enter the `postgres` password and type the SQL commands
     from [step 2 of the Linux installation](#postgresql).
   - MariaDB: installer from <https://mariadb.org/download/>. Then open **MySQL Client** (or
     **Command Prompt (MariaDB)** and `mariadb -u root -p`) and type the SQL commands from
     [step 2 of the Linux installation](#mariadb).

3. Download Bondamino and install the packages (choose the requirements file for your
   database; waitress is installed in addition):

   ```powershell
   git clone https://gitlab.com/d11m-group/bondamino.git C:\Bondamino
   cd C:\Bondamino
   py -m venv .venv
   .venv\Scripts\python.exe -m pip install -r requirements-postgresql.txt waitress
   ```

   For MariaDB use `requirements-mariadb.txt` (prebuilt `mysqlclient` packages exist for
   Windows), for SQLite `requirements.txt`.

4. Run the [guided setup](#the-guided-setup-in-detail):

   ```powershell
   cd C:\Bondamino\bondamino
   ..\.venv\Scripts\python.exe manage.py setup
   ```

   Answer as described for [Linux step 4](#step-4-run-the-guided-setup). For a quick test on the
   same computer add `localhost` to the domains, e.g. `bondamino.example.org,localhost`.

5. Test waitress:

   ```powershell
   cd C:\Bondamino\bondamino
   ..\.venv\Scripts\waitress-serve.exe --listen=127.0.0.1:8000 --trusted-proxy=127.0.0.1 --trusted-proxy-headers="x-forwarded-for x-forwarded-proto" bondamino.wsgi:application
   ```

   The window shows little or no output while waitress runs. Open <http://localhost:8000/login/>
   in a browser on the same computer; you should see the sign-in page (if you answered HTTPS with
   yes, signing in only works through the reverse proxy, see below). Stop with `Ctrl` + `C`.

   The two `--trusted-proxy` options are important behind a reverse proxy: by default waitress
   removes the `X-Forwarded-Proto` and `X-Forwarded-For` headers, so Bondamino would not know
   that the connection is HTTPS (invitation links would start with `http://`) and could not tell
   clients apart when limiting sign-in attempts. The options tell waitress to trust these headers
   from the proxy on the same machine.

6. Create a start script `C:\Bondamino\start-bondamino.cmd` (e.g. with `notepad
   C:\Bondamino\start-bondamino.cmd`) with this content. It writes the output, including errors
   and console emails, to `C:\Bondamino\logs\bondamino.log`:

   ```bat
   @echo off
   if not exist C:\Bondamino\logs mkdir C:\Bondamino\logs
   cd /d C:\Bondamino\bondamino
   C:\Bondamino\.venv\Scripts\waitress-serve.exe --listen=127.0.0.1:8000 --trusted-proxy=127.0.0.1 --trusted-proxy-headers="x-forwarded-for x-forwarded-proto" bondamino.wsgi:application >> C:\Bondamino\logs\bondamino.log 2>&1
   ```

7. Start Bondamino automatically when Windows starts. Choose one of the two options.

   **Option A: Task Scheduler** (built into Windows):

   ```powershell
   $action   = New-ScheduledTaskAction -Execute "C:\Bondamino\start-bondamino.cmd"
   $trigger  = New-ScheduledTaskTrigger -AtStartup
   $settings = New-ScheduledTaskSettingsSet -ExecutionTimeLimit ([TimeSpan]::Zero) -RestartCount 3 -RestartInterval (New-TimeSpan -Minutes 1) -AllowStartIfOnBatteries -DontStopIfGoingOnBatteries
   Register-ScheduledTask -TaskName "Bondamino" -Action $action -Trigger $trigger -Settings $settings -User "SYSTEM" -RunLevel Highest
   Start-ScheduledTask -TaskName "Bondamino"
   ```

   `-ExecutionTimeLimit ([TimeSpan]::Zero)` matters: otherwise Windows stops the task after three
   days. Stop and start it with `Stop-ScheduledTask -TaskName "Bondamino"` and
   `Start-ScheduledTask -TaskName "Bondamino"`.

   **Option B: a Windows service with NSSM.** Download NSSM from <https://nssm.cc/>, copy
   `nssm.exe` (from the `win64` folder) to `C:\Bondamino`, then:

   ```powershell
   C:\Bondamino\nssm.exe install Bondamino C:\Bondamino\start-bondamino.cmd
   C:\Bondamino\nssm.exe start Bondamino
   ```

   The service then appears in **Services** (`services.msc`) and starts with Windows. Stop and
   start it with `Stop-Service Bondamino` and `Start-Service Bondamino`.

   Check: <http://localhost:8000/login/> shows the sign-in page.

8. Put a reverse proxy with HTTPS in front. On Windows, [Caddy](#caddy) is the simplest choice;
   run it as a service with NSSM the same way (`nssm install Caddy C:\Caddy\caddy.exe run
   --config C:\Caddy\Caddyfile`). If your IT department runs IIS, it can act as reverse proxy with
   the **URL Rewrite** and **Application Request Routing** modules; it must preserve the original
   `Host` header and set `X-Forwarded-Proto: https`.

9. Open ports 80 and 443 in the Windows firewall for the reverse proxy (not port 8000):

   ```powershell
   New-NetFirewallRule -DisplayName "Bondamino HTTP/HTTPS" -Direction Inbound -Protocol TCP -LocalPort 80,443 -Action Allow
   ```

Continue with [After the installation](#after-the-installation).

**macOS** can run Bondamino standalone like Linux (gunicorn, with Python from python.org or
Homebrew), but it is rarely used as a server; for a Mac, the [Docker try-out](#try-out-with-docker)
is usually the better choice.

---

## Reverse proxy and HTTPS

Every installation that other people use should be reached via HTTPS through a reverse proxy.
Whatever proxy you use, it must:

1. pass on the original `Host` header (the domain in `ALLOWED_HOSTS`),
2. set `X-Forwarded-Proto` to `https`,
3. append the client address to `X-Forwarded-For` (Bondamino uses the **last** entry to limit
   sign-in attempts per client),
4. accept request bodies of at least 25 MB (uploads are limited to 20 MB; without a limit in the
   proxy, a too large file is only rejected after it has been received completely),
5. allow responses to take up to 120 seconds.

And in Bondamino set `HTTPS=1` (the setup question "served via HTTPS"). Bondamino itself does not
redirect HTTP to HTTPS and does not send HSTS; leave both to the proxy.

If there are two proxies in a row (e.g. a central load balancer of the institute in front of
nginx), the last `X-Forwarded-For` entry is the first proxy's address. Configure the inner proxy
to take the real client address from the outer one (in nginx: the `realip` module), otherwise
all clients share one address for the sign-in limits.

### nginx

See [step 7 of the Linux installation](#step-7-nginx-and-https).

### Apache

Enable the modules and create a site, e.g. `/etc/apache2/sites-available/bondamino.conf`:

```sh
sudo a2enmod proxy proxy_http headers ssl
```

```apache
<VirtualHost *:443>
    ServerName bondamino.example.org
    SSLEngine on
    SSLCertificateFile    /etc/ssl/certs/bondamino.crt
    SSLCertificateKeyFile /etc/ssl/private/bondamino.key

    ProxyPreserveHost On
    RequestHeader set X-Forwarded-Proto "https"
    ProxyPass        / http://127.0.0.1:8000/
    ProxyPassReverse / http://127.0.0.1:8000/
    ProxyTimeout 120
    LimitRequestBody 26214400
    Header always set Strict-Transport-Security "max-age=31536000"
</VirtualHost>

<VirtualHost *:80>
    ServerName bondamino.example.org
    Redirect permanent / https://bondamino.example.org/
</VirtualHost>
```

Apache adds `X-Forwarded-For` automatically. Enable the site with
`sudo a2ensite bondamino && sudo apachectl configtest && sudo systemctl reload apache2`. For a
Let's Encrypt certificate use `sudo apt install python3-certbot-apache` and
`sudo certbot --apache -d bondamino.example.org`.

### Caddy

Caddy (<https://caddyserver.com/>) gets and renews Let's Encrypt certificates automatically and
sets the forwarding headers itself. The whole configuration (`Caddyfile`):

```caddy
bondamino.example.org {
    request_body {
        max_size 25MB
    }
    header Strict-Transport-Security "max-age=31536000"
    reverse_proxy 127.0.0.1:8000
}
```

Automatic certificates need ports 80 and 443 to be reachable from the internet; inside an
institute network, add `tls /path/to/certificate.crt /path/to/private.key` with a certificate
from your IT department.

### Traefik and other proxies in Docker

Traefik and Caddy set the forwarding headers by default. Route the domain to the `bondamino`
service on port 8000, limit the request body (Traefik: the `buffering` middleware with
`maxRequestBodyBytes`), and set `HTTPS=1`.

---

## Docker Compose

Docker Compose runs Bondamino and, if needed, its database in containers on one host. The
repository contains one Compose file per database:

| File | Database | Containers |
|---|---|---|
| `docker-compose.yml` | PostgreSQL 17 (recommended) | `bondamino`, `db` |
| `docker-compose.mariadb.yml` | MariaDB 11 | `bondamino`, `db` |
| `docker-compose.sqlite.yml` | SQLite (file in the data volume) | `bondamino` |

On every start, the Bondamino container creates the secret key if needed (stored in the volume),
waits up to about 60 seconds for the database, applies database migrations, creates the cache
table, creates the first administrator if there is none yet, and starts gunicorn on port 8000.

### Step 1: Get the files

```sh
git clone https://gitlab.com/d11m-group/bondamino.git Bondamino
cd Bondamino
```

On a server, a folder such as `/opt/Bondamino` is a good place
(`sudo git clone https://gitlab.com/d11m-group/bondamino.git /opt/Bondamino`, then
`sudo chown -R $USER: /opt/Bondamino`).

### Step 2: Configure

Copy the template and edit it (`cp .env.example .env` and `nano .env` on Linux/macOS,
`Copy-Item .env.example .env` and `notepad .env` on Windows). Compose reads the file `.env` next
to the Compose file automatically. Write values without quotes.

| Variable | What to enter |
|---|---|
| `DB_PASSWORD` | Password for the database user (not needed for SQLite). **Letters and digits only**, because it becomes part of a URL. Generate one with `python3 -c "import secrets; print(secrets.token_hex(24))"`, `openssl rand -hex 24`, or in PowerShell `[guid]::NewGuid().ToString('N') + [guid]::NewGuid().ToString('N')`. It is set when the database is created on the first start; changing it later does not change the password inside the existing database (see [troubleshooting](#database-not-reachable)). |
| `ALLOWED_HOSTS` | The domain(s) under which Bondamino is reached, e.g. `bondamino.example.org`. The default `localhost,127.0.0.1` is enough for opening it on the host itself; add these two to your domain if you also want that on a server. |
| `HTTPS` | `1` when a reverse proxy provides HTTPS (normal on a server), `0` for plain `http://localhost`. |
| `ADMIN_USERNAME`, `ADMIN_EMAIL`, `ADMIN_PASSWORD` | The first administrator, created on the first start if there is no administrator yet. A weak password is accepted, with a warning in the log, so choose a long one. Remove `ADMIN_PASSWORD` after the first start (see [step 5](#step-5-after-the-first-start)). |
| `EMAIL_*`, `DEFAULT_FROM_EMAIL` | Your mail server (see the [reference](#email)). With `EMAIL_BACKEND` empty (the template's default), emails are only printed into the log. To send them, set `EMAIL_BACKEND=django.core.mail.backends.smtp.EmailBackend` and fill in `EMAIL_HOST`, `EMAIL_PORT`, `EMAIL_HOST_USER`, `EMAIL_HOST_PASSWORD` and `DEFAULT_FROM_EMAIL` (the template is prepared for STARTTLS on port 587). Note that `localhost` inside a container is the container itself, not the host. |
| `TIME_ZONE` | Default time zone for dates and times, e.g. `Europe/Berlin`; the default is `UTC`. Users can choose their own in their account settings. |
| `BONDAMINO_PORT` | Address and port on the host. The default `127.0.0.1:8000` makes Bondamino reachable only from the host itself, which is right for a reverse proxy on the same host and for a try-out. `8000` makes it reachable from other computers too (without a proxy); **never combine this with `HTTPS=1`**, because then clients could bypass the proxy and forge the `X-Forwarded-For` header used for limiting sign-in attempts. |
| `GUNICORN_WORKERS`, `GUNICORN_TIMEOUT` | Number of web server processes (default `3`) and request timeout in seconds (default `120`). |
| `LOG_LEVEL`, `CONN_MAX_AGE` | Amount of logging (default `WARNING`) and how long a database connection is reused (default `60` seconds); see the [reference](#configuration-reference). |
| `SECRET_KEY` | Leave empty: a key is generated on the first start and stored in the data volume. |

`DATABASE_URL` is set by the Compose file; the line in `.env` is ignored. `DEBUG` is passed and
must stay `0`. Only the variables listed in the `environment:` section of the Compose file reach
Bondamino; to set others (e.g. `HTTPS_PROXY`), see [Additional variables](#additional-variables).

### Step 3: Start

Build the image from the source code and start (only one of the three lines):

```sh
docker compose up -d --build                                  # PostgreSQL
docker compose -f docker-compose.mariadb.yml up -d --build    # MariaDB
docker compose -f docker-compose.sqlite.yml up -d --build     # SQLite
```

The first build takes a few minutes. A message that `registry.gitlab.com/d11m-group/bondamino` could not be
pulled is expected while the prebuilt image is not public; the image is built locally and tagged
with that name. If the maintainer has published the prebuilt image, you can leave out `--build`;
Compose then downloads the image instead of building it.

For MariaDB and SQLite, all further `docker compose` commands need the same `-f ...` option.

### Step 4: Check

```sh
docker compose ps
docker compose logs -f bondamino
```

`docker compose ps` should list `bondamino-bondamino-1` as `Up` (and `bondamino-db-1` as
`Up ... (healthy)`). The container names start with the project name, which is the folder name
(here `bondamino`). In the log, expect `Generated a new secret key in /data/secret_key` (first
start only), `Applying ...` lines (first start), `Superuser created successfully.` (first start)
and `Listening at: http://0.0.0.0:8000`. While the database starts up you may see
`Waiting for the database (1/30) ...`, which is normal. Stop following the log with `Ctrl` + `C`.

Open `http://localhost:8000/` on the host (if `localhost` is in `ALLOWED_HOSTS` and `HTTPS=0`),
or set up the [reverse proxy](#reverse-proxy-and-https) and open your domain. With `HTTPS=1`,
signing in over plain `http://` does not work by design.

### Step 5: After the first start

Remove the administrator password from the configuration, because it would otherwise stay
readable in the container settings: set `ADMIN_PASSWORD=` (empty) in `.env` and apply it with

```sh
docker compose up -d
```

(with `-f ...` for MariaDB/SQLite). Compose recreates the container with the new settings; the
data is kept. Continue with [After the installation](#after-the-installation).

### Everyday commands

| Task | Command |
|---|---|
| Show status | `docker compose ps` |
| Follow the log | `docker compose logs -f bondamino` |
| Stop (data is kept) | `docker compose down` |
| Start again | `docker compose up -d` |
| Restart | `docker compose restart bondamino` |
| Run a management command | `docker compose exec bondamino python manage.py <command>` |
| Create another administrator | `docker compose exec bondamino python manage.py createsuperuser` |
| Reset a password | `docker compose exec bondamino python manage.py changepassword <user name>` |
| Send a test email | `docker compose exec bondamino python manage.py sendtestemail you@example.org` |
| List the volumes | `docker volume ls` |

### Additional variables

The Compose files pass only selected variables to the container: those in `.env.example`
(including `DEBUG`, `LOG_LEVEL`, `CONN_MAX_AGE`, `GUNICORN_WORKERS` and `GUNICORN_TIMEOUT`, which
you simply set in `.env`). To set others, such as a proxy for outgoing connections
(`HTTPS_PROXY`), create a file `docker-compose.override.yml` next to `docker-compose.yml` instead
of editing the provided files (so that updates do not conflict):

```yaml
services:
  bondamino:
    environment:
      HTTPS_PROXY: http://proxy.example.org:3128
```

Compose reads `docker-compose.override.yml` automatically together with `docker-compose.yml`.
For the MariaDB and SQLite variants name both files explicitly, e.g.
`docker compose -f docker-compose.mariadb.yml -f docker-compose.override.yml up -d`.

### Storing the data in a host folder

By default the data lives in named Docker volumes (`<project>_bondamino-data` for uploads and the
secret key, `<project>_postgres-data` or `<project>_mariadb-data` for the database). To keep the
uploads in a folder of the host instead (e.g. for an existing backup system), override the
volume in `docker-compose.override.yml`:

```yaml
services:
  bondamino:
    volumes:
      - /srv/bondamino/data:/data
```

The folder must be writable for user ID 1000, the user Bondamino runs as in the container:
`sudo mkdir -p /srv/bondamino/data && sudo chown -R 1000:1000 /srv/bondamino/data`.

---

## Portainer

Portainer manages Docker through a web interface. Bondamino is deployed there as a **stack**
(Portainer's name for a Compose project). The environment in Portainer must be a Docker
standalone environment (not Swarm), because the image is built from the source code.

### Repository method (recommended)

Portainer downloads the repository and builds the image from the `Dockerfile`.

1. **Stacks → Add stack**, name: `bondamino` (lower case).
2. **Build method: Repository.**
   - Repository URL: `https://gitlab.com/d11m-group/bondamino`
   - Repository reference: `refs/heads/main`
   - Compose path: `docker-compose.yml` (PostgreSQL), `docker-compose.mariadb.yml` or
     `docker-compose.sqlite.yml`
3. **Environment variables:** add the variables from
   [Docker Compose step 2](#step-2-configure) (`DB_PASSWORD`, `ALLOWED_HOSTS`, `HTTPS`,
   `ADMIN_USERNAME`, `ADMIN_EMAIL`, `ADMIN_PASSWORD`, `EMAIL_*`, `DEFAULT_FROM_EMAIL`,
   `TIME_ZONE`, `BONDAMINO_PORT` and, if needed, `GUNICORN_WORKERS`, `GUNICORN_TIMEOUT`,
   `LOG_LEVEL`, `CONN_MAX_AGE`). **Advanced mode** lets you paste them as `NAME=value` lines;
   **Load variables from .env file** reads a prepared file. Also add

   ```ini
   BONDAMINO_IMAGE=bondamino-local:1
   ```

   This gives the self-built image a local name. Portainer builds an image only if it does not
   exist yet, so increasing this number is how you make it build a new one during an
   [update](#updating).
4. **Deploy the stack.** Building takes a few minutes. A message that the image could not be
   pulled is expected; the image is then built.
5. Check: **Containers → bondamino-bondamino-1 → Logs** should end with
   `Listening at: http://0.0.0.0:8000`. Bondamino listens on port 8000 of the host, by default
   only for the host itself (`127.0.0.1:8000`; see `BONDAMINO_PORT`).
6. Remove `ADMIN_PASSWORD` in the stack's environment variables and click **Update the stack**.

### Web editor method

The Web editor cannot build images, because it has no access to the source code. It therefore
needs an image that already exists:

- **Prebuilt image, if the maintainer has published it:** nothing to prepare; the Compose file
  uses `registry.gitlab.com/d11m-group/bondamino:latest`. If the registry is private, add the
  registry `registry.gitlab.com` with credentials (a GitLab deploy token or personal access
  token with the scope `read_registry`) under **Registries** first.
- **Self-built image:** build it once on the Docker host (SSH, Docker CLI):

  ```sh
  docker build -t bondamino-local:1 https://gitlab.com/d11m-group/bondamino.git#main
  ```

  and add `BONDAMINO_IMAGE=bondamino-local:1` to the environment variables below.

Then:

1. **Stacks → Add stack**, name `bondamino`, **Build method: Web editor**.
2. Paste the content of one of the Compose files (from the repository on GitLab) and **delete the
   line `build: .`**.
3. Add the environment variables as in the repository method, then **Deploy the stack**.

### Changing variables later

Open the stack, change the variables, and click **Update the stack** (Web editor) or
**Pull and redeploy** (Repository). The container is recreated; the data stays in the volumes.

---

## Other container platforms

The image runs on any platform that runs Linux containers (Kubernetes, Nomad, Synology Container
Manager, Unraid, cloud container services, ...). Build it and push it to a registry that your
platform can use:

```sh
docker build -t registry.example.org/bondamino:1 .
docker push registry.example.org/bondamino:1
```

(On a Mac with Apple silicon, build for typical servers with
`docker buildx build --platform linux/amd64 -t registry.example.org/bondamino:1 --push .`.)
If the maintainer has published the prebuilt image, `registry.gitlab.com/d11m-group/bondamino` (for
`linux/amd64` and `linux/arm64`) can be used directly.

| Property | Value |
|---|---|
| Port | `8000` (plain HTTP; put HTTPS in front) |
| Volume | `/data`: uploaded PDB files (`/data/media`), the generated secret key (`/data/secret_key`) and, without `DATABASE_URL`, the SQLite database (`/data/db.sqlite3`) |
| User | runs as user ID 1000 (`bondamino`); a mounted `/data` must be writable for it (Kubernetes: `fsGroup: 1000`) |
| Configuration | environment variables from the [reference](#configuration-reference); for PostgreSQL or MariaDB set `DATABASE_URL`, otherwise SQLite in `/data` is used |
| First administrator | `DJANGO_SUPERUSER_USERNAME`, `DJANGO_SUPERUSER_EMAIL`, `DJANGO_SUPERUSER_PASSWORD`; created only while there is no administrator. A weak password and a user name that already belongs to a normal user only produce a warning in the log; the container keeps running |
| Start-up | applies pending migrations on every start and creates the cache table; waits up to about 60 seconds for the database, then exits with "Database not reachable, giving up." |
| Other commands | arguments given to the container are run instead of gunicorn after the preparation, e.g. `python manage.py createsuperuser` |
| Web server tuning | `GUNICORN_WORKERS` (default `3`), `GUNICORN_TIMEOUT` (default `120` seconds) |

Example with plain `docker run` (SQLite in the volume):

```sh
docker volume create bondamino-data
docker run -d --name bondamino --restart unless-stopped \
  -p 127.0.0.1:8000:8000 -v bondamino-data:/data \
  -e ALLOWED_HOSTS=bondamino.example.org -e HTTPS=1 \
  -e DJANGO_SUPERUSER_USERNAME=admin -e DJANGO_SUPERUSER_EMAIL=you@example.org \
  -e DJANGO_SUPERUSER_PASSWORD=a-long-password \
  bondamino-local:1
```

**Running more than one replica** (several containers of the same installation) requires
PostgreSQL or MariaDB, a `/data` volume shared by all replicas (for the uploaded files) and the
same `SECRET_KEY` in all of them (set it explicitly instead of relying on the generated file).
Because every container applies migrations at start, start one replica first after an update.
For a working group, a single replica is sufficient.

---

## After the installation

### First sign-in

Open `https://bondamino.example.org/admin/` (or `http://localhost:8000/admin/` for a try-out) and
sign in with the administrator. This is the administration area (Django admin), where users,
work groups and languages are managed. The application itself is at
`https://bondamino.example.org/`. The admin follows the selected language; German labels are
given in brackets below.

### Creating users

There is no self-registration: users are created by an administrator and invited by email.

1. **Users (Benutzer) → Add.** Enter user name, email address (required), and optionally first
   and last name and work groups. **Save.**
2. Bondamino sends an invitation email with a link where the person sets their own password. The
   link is valid for **24 hours** and can be used once. The admin shows "The invitation was sent
   to ...", or an error message if sending failed (the user is created anyway).
3. The user list shows in the column **Password set (Passwort festgelegt)** who has accepted the
   invitation.
4. If the link has expired or the email did not arrive, select the user in the list, choose the
   action **Resend invitation (Einladung erneut senden)** and click **Go**. Alternatively, the
   person can use **Forgot password** on the sign-in page, which also works for invited users.

To give someone administrator rights, open the user and tick **Staff status** (access to the
admin) and **Superuser status** (all rights).

### Work groups

**Work groups (Arbeitsgruppen) → Add:** a name, a leader and the members. The leader is always
added as a member. Files and notes can be shared with work groups; users can also be added to
groups on their user page. How groups are used is described in the [tutorial](TUTORIAL.md).

### Languages

German and English are created during installation and are offered in the language selection.
Under **Languages (Sprachen)** you can rename or remove them (without any language, the language
field disappears from the users' account settings and the default language applies). Only German and English
translations are included; other languages would need translation files (see
[DOCUMENTATION.md](DOCUMENTATION.md)).

### Checking email delivery

- Send a test email with the current configuration:

  ```sh
  python manage.py sendtestemail you@example.org                              # standalone, in the folder bondamino
  docker compose exec bondamino python manage.py sendtestemail you@example.org  # Docker
  ```

  An error message (e.g. `SMTPAuthenticationError`, `Connection refused`) shows what is wrong.
- Then create a test user with your own email address and check that the invitation arrives
  (also in the spam folder).
- With the console backend (`EMAIL_BACKEND` empty, the default, or `...console.EmailBackend`), emails are
  written to the log instead: `sudo journalctl -u bondamino` (systemd),
  `C:\Bondamino\logs\bondamino.log` (Windows), `docker compose logs bondamino` (Docker). For a
  try-out you can copy the invitation link from there.

### Checking HTTPS and the proxy

- `HTTPS=1` must be set when users reach Bondamino through HTTPS. Then cookies are only sent
  over HTTPS, and the proxy's `X-Forwarded-Proto` and `X-Forwarded-For` headers are trusted.
- Invitation links must start with `https://` and your domain. If they start with `http://` or
  contain `127.0.0.1`, the proxy does not pass the `Host` or `X-Forwarded-Proto` header (or, on
  Windows, waitress lacks the `--trusted-proxy` options).
- Make sure Bondamino's port 8000 is not reachable from outside (bound to `127.0.0.1`, the
  default of `BONDAMINO_PORT` in Docker Compose), so that nobody can bypass the proxy and forge
  the `X-Forwarded-For` header.
- Standalone: `python manage.py check --deploy` lists security settings. The warnings
  `security.W004` (HSTS) and `security.W008` (redirect to HTTPS) are expected, because the proxy
  handles both. Other warnings should be looked at.

### Removing the initial administrator password

With Docker, Portainer and other container platforms, remove `ADMIN_PASSWORD`
(`DJANGO_SUPERUSER_PASSWORD`) from the configuration after the first start, as described in
[Docker Compose step 5](#step-5-after-the-first-start). Changing it later has no effect on the
existing administrator; use `changepassword` or the admin to change passwords.

---

## Configuration reference

Bondamino reads its configuration from environment variables. In a standalone installation, the
file `bondamino/.env` (next to `manage.py`) is read in addition; variables that are already set
in the environment take precedence over the file. With Docker Compose, the `.env` file next to
the Compose file provides values for the Compose file, which passes them to the container.

Format rules: one `NAME=value` per line, lines starting with `#` are comments, values are written
without quotes. An empty value counts as not set, so the default applies. Yes/no values accept
`1`, `true`, `yes` or `on` (anything else means no). Lists are comma-separated.

### General

| Variable | Default | Meaning |
|---|---|---|
| `SECRET_KEY` | none (Docker: generated) | Long random value that signs sessions and password links. Keep it secret. Bondamino does not start without it. The setup generates one; in Docker it is generated on the first start and stored in `/data/secret_key`. Changing it signs out all users and invalidates open invitation links, but not the passwords. |
| `SECRET_KEY_FILE` | none (Docker: `/data/secret_key`) | File that contains the secret key; used if `SECRET_KEY` is empty. |
| `DEBUG` | `0` | `1` = development mode with detailed error pages. Never on a server that others can reach. |
| `ALLOWED_HOSTS` | with `DEBUG=1`: `localhost,127.0.0.1,[::1]`; otherwise none (Compose files: `localhost,127.0.0.1`) | Host names under which Bondamino is reached, comma-separated, e.g. `bondamino.example.org,localhost`. A leading dot allows all subdomains (`.example.org`). Requests for other names get "400 Bad Request". Links in emails use the name under which the page was requested, so the public domain must be listed and the proxy must pass it on. |
| `HTTPS` | `0` | `1` = Bondamino is served via HTTPS through a reverse proxy: session and CSRF cookies are only sent over HTTPS, `X-Forwarded-Proto: https` marks a request as secure, and the client address for limiting sign-in attempts is taken from the last entry of `X-Forwarded-For`. |
| `CSRF_TRUSTED_ORIGINS` | derived from `ALLOWED_HOSTS` (`https://<host>` with `HTTPS=1`, otherwise `http://<host>`) | Origins from which forms may be submitted. Only needed if the proxy changes the host name or port, e.g. `https://bondamino.example.org`. Must include the scheme. |
| `TIME_ZONE` | `UTC` | Default time zone for displaying dates and times, as an IANA name, e.g. `Europe/Berlin`. Each user can choose their own time zone in their account settings; this one applies to everybody who has not. |
| `LOG_LEVEL` | `WARNING` | Level of the console log: `DEBUG`, `INFO`, `WARNING`, `ERROR`. |
| `BONDAMINO_ENV_FILE` | `bondamino/.env` | Standalone: path of the configuration file to read instead. |

### Database

| Variable | Default | Meaning |
|---|---|---|
| `DATABASE_URL` | `sqlite:///db.sqlite3` (Docker image: `sqlite:////data/db.sqlite3`) | Database connection, see below. Set by the Compose files. |
| `CONN_MAX_AGE` | `60` | Seconds a database connection is kept open and reused; `0` closes it after every request. |

#### DATABASE_URL

| Database | Format | Example |
|---|---|---|
| SQLite, relative path | `sqlite:///<file>` (relative to the folder of `manage.py`) | `sqlite:///db.sqlite3` |
| SQLite, absolute path (Linux/macOS) | `sqlite:////<path>` (four slashes) | `sqlite:////var/lib/bondamino/db.sqlite3` |
| SQLite, absolute path (Windows) | `sqlite:///<drive>:/<path>` | `sqlite:///C:/Bondamino-data/db.sqlite3` |
| PostgreSQL | `postgres://<user>:<password>@<host>:<port>/<database>` (also `postgresql://`) | `postgres://bondamino:secret@localhost:5432/bondamino` |
| MariaDB / MySQL | `mariadb://<user>:<password>@<host>:<port>/<database>` (also `mysql://`) | `mariadb://bondamino:secret@localhost:3306/bondamino` |

Options for the driver can be added after `?`, e.g.
`postgres://bondamino:secret@db.example.org:5432/bondamino?sslmode=require`. For MariaDB/MySQL
the character set is `utf8mb4` unless `?charset=...` is given. The port can be left out to use
the default.

**Special characters in user name or password** must be URL-encoded, because characters such as
`@`, `:` and `/` have a meaning in a URL. The guided setup does this automatically. By hand:

| Character | Encoded | Character | Encoded |
|---|---|---|---|
| `@` | `%40` | `#` | `%23` |
| `:` | `%3A` | `?` | `%3F` |
| `/` | `%2F` | `&` | `%26` |
| `%` | `%25` | space | `%20` |

Example: the password `p@ss:word` becomes `p%40ss%3Aword`. Python encodes a password for you:

```sh
python3 -c "from urllib.parse import quote; print(quote(input('Password: '), safe=''))"
```

The Compose files insert `DB_PASSWORD` into the URL without encoding; that is why it must consist
of letters and digits only.

### Files

| Variable | Default | Meaning |
|---|---|---|
| `MEDIA_ROOT` | `bondamino/media` (Docker image: `/data/media`) | Folder for uploaded PDB files. They are never served directly, only through Bondamino after checking the access rights. Include it in backups. |
| `STATIC_ROOT` | `bondamino/staticfiles` (Docker image: `/app/staticfiles`) | Target of `manage.py collectstatic`. Bondamino serves these files itself (WhiteNoise); the reverse proxy does not need to. |

### Email

Emails are used for invitations and forgotten passwords.

| Variable | Default | Meaning |
|---|---|---|
| `EMAIL_BACKEND` | `django.core.mail.backends.console.EmailBackend` (also when empty, as in `.env.example`) | `django.core.mail.backends.smtp.EmailBackend` sends emails via SMTP with the values below; the console backend only prints them into the log. |
| `EMAIL_HOST` | `localhost` | SMTP server |
| `EMAIL_PORT` | `25` | SMTP port: usually `587` (STARTTLS), `465` (SSL/TLS) or `25` |
| `EMAIL_HOST_USER` | empty | SMTP user name; empty if the server needs no login |
| `EMAIL_HOST_PASSWORD` | empty | SMTP password |
| `EMAIL_USE_TLS` | `0` | `1` = STARTTLS (usually port 587) |
| `EMAIL_USE_SSL` | `0` | `1` = SSL/TLS (usually port 465). Do not set together with `EMAIL_USE_TLS`. |
| `DEFAULT_FROM_EMAIL` | `Bondamino <noreply@localhost>` | Sender of all emails, e.g. `Bondamino <bondamino@example.org>`. Use an address your mail server accepts as sender. |

### Docker and Compose only

| Variable | Default | Meaning |
|---|---|---|
| `DB_PASSWORD` | none (required for PostgreSQL/MariaDB) | Compose: password of the database user; letters and digits only. Only used when the database volume is created. |
| `BONDAMINO_PORT` | `127.0.0.1:8000` | Compose: address and port on the host. The default is only reachable from the host itself (for a reverse proxy on the same host). `8000` makes Bondamino reachable from other computers without a proxy; never use that together with `HTTPS=1`, because clients could then forge the `X-Forwarded-For` header. |
| `BONDAMINO_IMAGE` | `registry.gitlab.com/d11m-group/bondamino:latest` | Compose: name of the image to use or to build, e.g. `bondamino-local:1`. |
| `ADMIN_USERNAME`, `ADMIN_EMAIL`, `ADMIN_PASSWORD` | empty | Compose: first administrator; passed to the container as the three variables below. |
| `DJANGO_SUPERUSER_USERNAME`, `DJANGO_SUPERUSER_EMAIL`, `DJANGO_SUPERUSER_PASSWORD` | empty | Container: first administrator, created at start if user name and password are set and no administrator exists yet. A weak password is accepted with a warning in the log. If the user name already belongs to a normal user, nothing is created and a warning is logged. |
| `GUNICORN_WORKERS` | `3` | Container: number of web server processes. |
| `GUNICORN_TIMEOUT` | `120` | Container: seconds before a request is aborted. |

The Compose files pass these variables to the container: `ALLOWED_HOSTS`, `CSRF_TRUSTED_ORIGINS`,
`HTTPS`, `TIME_ZONE`, `DEBUG`, `LOG_LEVEL`, `CONN_MAX_AGE`, `GUNICORN_WORKERS`,
`GUNICORN_TIMEOUT`, `SECRET_KEY`, all `EMAIL_*` variables, `DEFAULT_FROM_EMAIL` and the
administrator variables. `DATABASE_URL` is fixed in each Compose file.

---

## Updating

Read the release notes or the commit history first, and **make a [backup](#backup-and-restore)
before every update**. Database migrations are applied automatically or with `migrate`; they
cannot be undone without the backup.

### Standalone (Linux)

```sh
sudo -u bondamino -H bash
cd /opt/Bondamino
git pull
. .venv/bin/activate
pip install -r requirements-postgresql.txt     # the file of your database
cd bondamino
python manage.py migrate
python manage.py createcachetable
python manage.py collectstatic --noinput
exit
sudo systemctl restart bondamino
```

What the commands do: `git pull` fetches the new version (expected: a list of changed files, or
`Already up to date.`), `pip install` installs new or updated packages, `migrate` updates the
database tables, `createcachetable` creates the cache table if a new version needs it (no change
if it exists), and `collectstatic` updates styles and scripts. Instead of the last three
commands you can run `python manage.py setup --finish`.

Check with `sudo systemctl status bondamino` and by opening the site; if styles are missing,
reload the page with `Ctrl` + `F5`.

### Standalone (Windows)

Stop Bondamino (`Stop-ScheduledTask -TaskName "Bondamino"` or `Stop-Service Bondamino`), then:

```powershell
cd C:\Bondamino
git pull
.venv\Scripts\python.exe -m pip install -r requirements-postgresql.txt waitress
cd bondamino
..\.venv\Scripts\python.exe manage.py setup --finish
```

and start it again.

### Installed from a ZIP file

Download the new ZIP, extract it, and copy over the old folder everything except your data:
keep `bondamino/.env`, `bondamino/media`, `bondamino/db.sqlite3` (SQLite) and `.venv`. Then
continue with `pip install` and the following commands above.

### Docker Compose

```sh
cd /opt/Bondamino
git pull
docker compose up -d --build
```

(with `-f docker-compose.mariadb.yml` or `-f docker-compose.sqlite.yml` for the other variants).
Compose rebuilds the image with the new source code and recreates the container; the container
applies the migrations at start. Check with `docker compose logs -f bondamino`.

Also update the base system regularly, even without a new Bondamino version:

```sh
docker compose build --pull     # newest Python base image with security fixes
docker compose pull db          # newest minor version of the database image
docker compose up -d
docker image prune              # remove old, unused images
```

Do not change the major version of the database image (e.g. `postgres:17` to `postgres:18`) by
editing the Compose file: the existing data must then be moved with a dump and restore.

If the prebuilt image is published and you use it instead of building, update with
`docker compose pull` followed by `docker compose up -d`.

### Portainer

- **Repository method:** open the stack, increase `BONDAMINO_IMAGE` (e.g. from
  `bondamino-local:1` to `bondamino-local:2`) and click **Pull and redeploy**. Portainer fetches
  the new source code and, because an image with the new name does not exist yet, builds it.
  Remove old images under **Images** (filter "unused").
- **Web editor with a self-built image:** build a new image on the host
  (`docker build -t bondamino-local:2 https://gitlab.com/d11m-group/bondamino.git#main`), change
  `BONDAMINO_IMAGE` in the stack and click **Update the stack**.
- **Web editor with the prebuilt image:** click **Update the stack** with **Re-pull image**
  enabled.

---

## Backup and restore

Back up three things:

1. **the database** (users, groups, notes, file metadata),
2. **the uploaded PDB files** (`MEDIA_ROOT`; in Docker the volume `/data`, which also contains the
   secret key and, for SQLite, the database),
3. **the configuration**: `bondamino/.env` (standalone) or the Compose `.env` file / the Portainer
   variables. It contains passwords and the secret key: store it securely.

Database and files belong together; back them up at the same time and keep several generations
(e.g. daily for a week, weekly for a month), preferably on another machine. Test a restore from
time to time.

### Standalone

PostgreSQL (asks for the password of the database user):

```sh
pg_dump -h localhost -U bondamino bondamino > bondamino-db.sql
```

MariaDB:

```sh
mariadb-dump --single-transaction -u bondamino -p bondamino > bondamino-db.sql
```

SQLite: stop Bondamino (`sudo systemctl stop bondamino`) and copy
`/opt/Bondamino/bondamino/db.sqlite3`, then start it again.

Uploaded files and configuration:

```sh
sudo tar czf bondamino-files.tgz -C /opt/Bondamino/bondamino media .env
```

On Windows, stop Bondamino and copy the folder `C:\Bondamino\bondamino\media`, the file
`C:\Bondamino\bondamino\.env` and, for SQLite, `db.sqlite3`; for PostgreSQL/MariaDB run the
dump commands above in PowerShell (they are installed with the database server).

### Docker Compose

Find the volume names with `docker volume ls`; they start with the project name (the folder
name, e.g. `bondamino_bondamino-data`).

PostgreSQL:

```sh
docker compose exec -T db pg_dump -U bondamino bondamino > bondamino-db.sql
```

MariaDB:

```sh
docker compose -f docker-compose.mariadb.yml exec -T db sh -c 'mariadb-dump --single-transaction -u bondamino -p"$MARIADB_PASSWORD" bondamino' > bondamino-db.sql
```

(`-T` prevents line-ending problems in the output file.)

Data volume (uploads, secret key, and for SQLite the database; for SQLite stop the container
first with `docker compose -f docker-compose.sqlite.yml stop bondamino`):

```sh
docker run --rm -v bondamino_bondamino-data:/data -v "$PWD":/backup alpine tar czf /backup/bondamino-data.tgz -C /data .
```

(In PowerShell, write `"${PWD}:/backup"` instead of `"$PWD":/backup`.) Also keep a copy of the
Compose `.env` file.

### Restore

Restore into an installation of the **same Bondamino version** as the backup; update afterwards.

Standalone:

1. Stop Bondamino: `sudo systemctl stop bondamino`.
2. Empty the database and import the dump.

   PostgreSQL:

   ```sh
   sudo -u postgres dropdb bondamino
   sudo -u postgres createdb -O bondamino bondamino
   psql -h localhost -U bondamino -d bondamino < bondamino-db.sql
   ```

   MariaDB:

   ```sh
   sudo mariadb -e 'DROP DATABASE bondamino; CREATE DATABASE bondamino CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;'
   mariadb -u bondamino -p bondamino < bondamino-db.sql
   ```

   SQLite: copy the backed-up `db.sqlite3` back.
3. Restore the files and configuration and fix the owner:

   ```sh
   sudo tar xzf bondamino-files.tgz -C /opt/Bondamino/bondamino
   sudo chown -R bondamino:bondamino /opt/Bondamino/bondamino
   ```

4. Start Bondamino: `sudo systemctl start bondamino`.

Docker Compose (PostgreSQL; replace `bondamino_` by your project name):

```sh
docker compose stop bondamino
docker compose exec -T db dropdb -U bondamino bondamino
docker compose exec -T db createdb -U bondamino bondamino
docker compose exec -T db psql -U bondamino -d bondamino < bondamino-db.sql
docker run --rm -v bondamino_bondamino-data:/data -v "$PWD":/backup alpine sh -c "rm -rf /data/* && tar xzf /backup/bondamino-data.tgz -C /data"
docker compose start bondamino
```

MariaDB: as above with `-f docker-compose.mariadb.yml`, but instead of `dropdb`, `createdb` and
`psql`:

```sh
docker compose -f docker-compose.mariadb.yml exec -T db sh -c 'mariadb -u bondamino -p"$MARIADB_PASSWORD" bondamino' < bondamino-db.sql
```

(the dump replaces the existing tables). SQLite: only restore the data volume.

On a new host, start the stack once with the same `DB_PASSWORD` (so that the database and volumes
exist), then restore as above.

---

## Moving to another server

1. Install Bondamino on the new server with the same method (or another one) and the **same
   version**, but do not create users there.
2. On the old server, stop Bondamino and make a [backup](#backup-and-restore).
3. Copy the backup files to the new server (e.g. `scp bondamino-* user@new-server:`).
4. Restore them on the new server. Keep the old `SECRET_KEY` (in `.env` or in the data volume)
   if users should stay signed in and open invitation links should stay valid; otherwise a new
   key is fine.
5. Adjust `ALLOWED_HOSTS` if the domain changes, move the DNS entry to the new server and set up
   the reverse proxy and certificate.
6. Check sign-in, a few files and email delivery, then shut down the old installation.

**Changing the database type** (e.g. from SQLite to PostgreSQL): a database dump cannot be
imported into another database type. Use Django's data export instead, and test it on a copy
first:

```sh
# in the old installation (standalone; with Docker prefix "docker compose exec bondamino")
python manage.py dumpdata --natural-foreign --natural-primary -e contenttypes -e auth.permission -e admin.logentry -e sessions --indent 2 -o bondamino-data.json

# in the new installation with the new, empty database (tables created with migrate or the setup)
python manage.py loaddata bondamino-data.json
```

Copy the uploaded files (`MEDIA_ROOT`) to the new installation as well.

---

## Uninstalling

Make a [backup](#backup-and-restore) first if you might need the data again.

**Standalone (Linux):**

```sh
sudo systemctl disable --now bondamino
sudo rm /etc/systemd/system/bondamino.service
sudo systemctl daemon-reload
sudo rm /etc/nginx/sites-enabled/bondamino /etc/nginx/sites-available/bondamino
sudo systemctl reload nginx
sudo rm -rf /opt/Bondamino                      # code, configuration, uploads, SQLite database
sudo userdel -r bondamino
sudo -u postgres dropdb bondamino               # PostgreSQL
sudo -u postgres dropuser bondamino
# MariaDB: sudo mariadb -e "DROP DATABASE bondamino; DROP USER 'bondamino'@'localhost';"
```

To remove the certificate as well: `sudo certbot delete --cert-name bondamino.example.org`.

**Standalone (Windows):** stop and remove the task (`Unregister-ScheduledTask -TaskName
"Bondamino"`) or service (`C:\Bondamino\nssm.exe remove Bondamino confirm`), delete
`C:\Bondamino`, remove the firewall rule
(`Remove-NetFirewallRule -DisplayName "Bondamino HTTP/HTTPS"`) and drop the database as above.

**Docker Compose:**

```sh
docker compose down          # removes the containers, keeps the data volumes
docker compose down -v       # removes the containers AND all data (database, uploads)
docker image rm registry.gitlab.com/d11m-group/bondamino:latest
```

(with `-f ...` for MariaDB/SQLite; use the image name from `BONDAMINO_IMAGE` if you changed it.)
Then delete the folder with the repository.

**Portainer:** **Stacks →** select the stack **→ Remove**. Remove the volumes under **Volumes**
(`bondamino_bondamino-data` and the database volume) and the image under **Images**.

---

## Security notes for operators

- **Sign-in limits.** After **5 failed attempts for one user name from one client address**, or
  **20 failed attempts from one client address in total**, sign-in from that address is blocked
  for **15 minutes**. This also applies to `/admin/`. While blocked, even the correct password is
  rejected. A successful sign-in resets only the counter for that user name. The forgotten-password
  form sends an email for at most **3 requests per email address** and **10 requests per client
  address** within an hour; further requests get the same answer but no email. The counters are stored in the
  database cache table (`createcachetable`); if the table is missing, nothing is limited and an
  error is logged. Behind a proxy, set `HTTPS=1` and let the proxy set `X-Forwarded-For`,
  otherwise all users share the proxy's address and one person's failed attempts can block
  everybody.
- **Keep the secret key secret.** Anyone who knows `SECRET_KEY` can forge sessions and password
  links. Do not commit `.env` to Git or send it by email; it is excluded in `.gitignore`. If it
  leaks, generate a new key (`python -c "import secrets; print(secrets.token_urlsafe(50))"`) and
  restart; all users are signed out.
- **DEBUG off.** `DEBUG=1` shows internal information on error pages. Never use it on a server.
- **HTTPS only.** Run Bondamino behind a reverse proxy with HTTPS, redirect HTTP to HTTPS and send
  HSTS from the proxy. Bind Bondamino's port 8000 to `127.0.0.1` or block it in the firewall.
  With Docker, published ports bypass host firewalls such as `ufw`, so keep the default
  `BONDAMINO_PORT=127.0.0.1:8000`. If you set `BONDAMINO_PORT=8000` to reach Bondamino from other
  computers without a proxy, keep `HTTPS=0`: with `HTTPS=1` clients could forge the
  `X-Forwarded-For` header and get around the sign-in limits.
- **Upload size.** Limit the request body in the proxy to about 25 MB (Bondamino accepts PDB files
  up to 20 MB).
- **Passwords and invitations.** Invitation and reset links are valid for 24 hours and only once.
  Remove `ADMIN_PASSWORD` after the first start. Give administrator rights only to the people who
  manage users.
- **Uploaded files** are only delivered to users who may see them; never configure the proxy to
  serve the media folder directly.
- **No external services.** All scripts, styles and fonts are delivered by Bondamino itself. The
  only outgoing connection is to `https://files.rcsb.org` when a structure is loaded by PDB ID.
  See [TRANSPARENCY.md](TRANSPARENCY.md).
- **Keep everything updated:** the operating system (`sudo apt update && sudo apt upgrade`,
  or automatic security updates with `unattended-upgrades`), the database server, the reverse
  proxy, Bondamino itself and, with Docker, the images (see [Updating](#updating)).
- **Backups** are part of security: keep them on another machine and protect them like the
  server, because they contain all data and the configuration.
- **Database access.** The database should only be reachable from Bondamino (PostgreSQL and
  MariaDB listen only on `localhost` by default; the Compose files do not publish the database
  port).

---

## Troubleshooting

First look at the log; it usually names the problem:

| Installation | Log |
|---|---|
| Standalone Linux | `sudo journalctl -u bondamino -n 100` |
| Standalone Windows | `C:\Bondamino\logs\bondamino.log` |
| Docker Compose | `docker compose logs --tail 100 bondamino` |
| Portainer | **Containers →** `...-bondamino-1` **→ Logs** |
| Reverse proxy | nginx: `/var/log/nginx/error.log`; Apache: `/var/log/apache2/error.log` |

For more detail, set `LOG_LEVEL=INFO` and restart.

### "Bad Request (400)"

The address in the browser is not in `ALLOWED_HOSTS`. Add the exact host name (without
`https://` and without port), e.g. `ALLOWED_HOSTS=bondamino.example.org`, and restart. Common
cases: opening `http://127.0.0.1:8000` when only `localhost` is listed (or vice versa); opening
the server by its IP address; a reverse proxy that does not pass the original `Host` header
(nginx: `proxy_set_header Host $host;`, Apache: `ProxyPreserveHost On`).

### "Forbidden (403) CSRF verification failed"

- The site is opened via HTTPS but `HTTPS=0`: set `HTTPS=1` and make sure the proxy sends
  `X-Forwarded-Proto: https`.
- `HTTPS=1` but the site is opened via plain `http://` (e.g. directly on port 8000): the cookies
  are then not stored ("CSRF cookie not set"). Use the HTTPS address.
- The proxy changes the host name or port: set `CSRF_TRUSTED_ORIGINS` to the address in the
  browser, e.g. `https://bondamino.example.org` or `https://bondamino.example.org:8443`.
- Cookies are blocked in the browser, or the form was open for a long time: reload and retry.

### Emails do not arrive

1. Check the configuration with `sendtestemail` (see
   [Checking email delivery](#checking-email-delivery)). If nothing is sent and the email
   appears in the log, the console backend is active: set
   `EMAIL_BACKEND=django.core.mail.backends.smtp.EmailBackend`.
2. `Connection refused` / timeout: wrong `EMAIL_HOST` or `EMAIL_PORT`, or a firewall blocks
   outgoing mail. In Docker, `localhost` is the container itself; use the host name of the
   mail server.
3. `SSL` / `wrong version number`: port and encryption do not match (587 with `EMAIL_USE_TLS=1`,
   465 with `EMAIL_USE_SSL=1`; never both).
4. `Authentication failed`: check `EMAIL_HOST_USER` and `EMAIL_HOST_PASSWORD`; some providers
   require an app password.
5. `Sender address rejected`: set `DEFAULT_FROM_EMAIL` to an address the mail server may send
   from.
6. Sent but not received: check the spam folder, and ask your IT department whether the sender
   domain needs SPF/DKIM entries.

Users whose invitation did not arrive can get a new one with the admin action **Resend
invitation**.

### Database not reachable

- Standalone: `python manage.py setup` tests the connection and shows the error. Check that the
  database server runs (`sudo systemctl status postgresql` or `mariadb`), and host, port, name,
  user and password in `DATABASE_URL`. Special characters in the password must be
  [URL-encoded](#database_url).
- `password authentication failed` (PostgreSQL) or `Access denied` (MariaDB) in Docker after
  changing `DB_PASSWORD`: the password is only set when the database volume is created. Either
  restore the old `DB_PASSWORD`, or change the password inside the database, e.g.
  `docker compose exec db psql -U bondamino -d bondamino -c "ALTER USER bondamino PASSWORD 'new-password';"`.
  For a fresh test installation without data, `docker compose down -v` deletes the volumes so
  that they are created again (this deletes all data).
- Docker log `Database not reachable, giving up.`: the database container did not start; check
  `docker compose logs db`.
- `No module named 'psycopg'` / `'MySQLdb'` / "Error loading ... module": the driver is not
  installed; run `pip install -r requirements-postgresql.txt` or `requirements-mariadb.txt` in
  the virtual environment.

### Missing cache table

Log message `Throttle cache not available; run "manage.py createcachetable"` or an error about
`bondamino_cache`: the cache table for the sign-in limits is missing (usually after an update
without `createcachetable`). Run `python manage.py createcachetable` (standalone, in the folder
`bondamino`) and restart. Docker creates it automatically on every start.

### Port already in use

`Address already in use` / `Connection in use: ('127.0.0.1', 8000)` / "port is already
allocated": another program uses port 8000, often a second Bondamino or another development
server. Find it with `sudo ss -tlnp | grep 8000` (Linux), `lsof -i :8000` (macOS) or
`Get-NetTCPConnection -LocalPort 8000` (Windows), stop it, or use another port: standalone
`--bind 127.0.0.1:8001` (gunicorn), `--listen=127.0.0.1:8001` (waitress) or
`python manage.py runserver 8001`, and adjust the proxy; Docker `BONDAMINO_PORT=127.0.0.1:8001`.

### Permission errors

- Docker: `Permission denied` for `/data/...` (e.g. when writing `/data/secret_key` or
  uploading): the mounted host folder is not writable for user ID 1000. Run
  `sudo chown -R 1000:1000 /path/to/folder`. On NAS systems, give user ID 1000 write access in
  the folder settings.
- Standalone: the service user cannot read `.env` or write `media`, e.g. because the setup ran as
  another user. Run `sudo chown -R bondamino:bondamino /opt/Bondamino`. If Bondamino does not
  start and the log says "The SECRET_KEY setting must not be empty", `.env` was not read: check
  its location (next to `manage.py`) and its owner.

### Page without styles (CSS missing)

The static files were not collected. Standalone: run `python manage.py collectstatic --noinput`
in the folder `bondamino` and restart. It is also run by the setup and should be part of every
update. In Docker they are collected when the image is built; rebuild with `--build`. Finally
reload the page with `Ctrl` + `F5`.

### mysqlclient cannot be installed

`pip install -r requirements-mariadb.txt` fails with errors such as `Can not find valid
pkg-config name`, `mysql_config not found` or `Python.h: No such file or directory`: the
compiler or the client library headers are missing.

- Debian/Ubuntu: `sudo apt install build-essential pkg-config python3-dev default-libmysqlclient-dev`
- Fedora/RHEL: `sudo dnf install gcc pkgconf-pkg-config python3-devel mariadb-connector-c-devel`
- macOS (Homebrew): `brew install pkg-config mysql-client`, then
  `export PKG_CONFIG_PATH="$(brew --prefix mysql-client)/lib/pkgconfig"` before `pip install`
- Windows: prebuilt packages exist for current Python versions; update pip
  (`python -m pip install --upgrade pip`) and try again.

Then repeat the `pip install`. Alternatively use PostgreSQL, whose driver needs no compiler, or
Docker, whose image contains both drivers.

### Other problems

- **Loading a structure by PDB ID fails** on a server: the server needs outgoing HTTPS access to
  `files.rcsb.org`. If your institute uses a web proxy, set `HTTPS_PROXY=http://proxy:port` for
  Bondamino (standalone: as a line in `bondamino/.env`; Docker: see
  [Additional variables](#additional-variables)) and restart.
- **Upload fails with "413 Request Entity Too Large"**: raise the request size limit in the
  proxy (`client_max_body_size 25m;` in nginx).
- **"502 Bad Gateway" or "504 Gateway Timeout"** from the proxy: Bondamino is not running or not
  listening on the address the proxy uses; check the service and the log.
- **Docker container restarts repeatedly**: look at `docker compose logs bondamino`; typical
  causes are an unreachable database or missing write permission on `/data`.
- **No administrator was created (Docker)**: the log shows a warning, e.g. that the user
  `ADMIN_USERNAME` already exists but is no administrator. Give that user administrator rights
  (`docker compose exec bondamino python manage.py createsuperuser` creates another one), or
  choose another `ADMIN_USERNAME` and restart with `docker compose up -d`.
- **Not reachable from other computers (Docker)**: by default the port is bound to `127.0.0.1`,
  i.e. only the host itself can open Bondamino. Put a [reverse proxy](#reverse-proxy-and-https)
  in front, or set `BONDAMINO_PORT=8000` (only with `HTTPS=0`) and add the server's name or IP
  address to `ALLOWED_HOSTS`.
- **Signed out after every restart (Docker)**: the secret key is not kept because `/data` is not
  a persistent volume. Mount a volume or set `SECRET_KEY`.
