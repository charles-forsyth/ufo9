## Mounting HPCC Filesystems on Your Local Machine with SSHFS

**Introduction:**

This article explains how to securely mount a directory from the UCR HPCC (High-Performance Computing Cluster) onto your local laptop using SSHFS (SSH Filesystem). This allows you to access and modify files on the cluster as if they were stored locally, providing a convenient way to transfer and manage data.

**Table of Contents:**

1.  Prerequisites
2.  Installing SSHFS
3.  Establishing the SSHFS Connection
4.  Using the Mounted Filesystem
5.  Disconnecting the SSHFS Connection
6.  Troubleshooting
7.  Conclusion & Additional Resources

**1. Prerequisites:**

Before you begin, ensure you have the following:

*   **HPCC Account:**  You must have an active account on the HPCC cluster.
*   **Web Consol SSH Access:**  You can access the HPCC head nodes through our web console. This is the preferred method of SSH access.
*   **SSHFS Installed on your Local Machine:** SSHFS needs to be installed on your local machine (laptop).
*   **Local Directory:**  Create a local directory on your laptop where you want to mount the HPCC filesystem.

**2. Installing SSHFS:**

The installation process varies depending on your operating system. Here are instructions for some common operating systems:

*   **Linux (Debian/Ubuntu):**

    ```bash
    sudo apt-get update
    sudo apt-get install sshfs
    ```

*   **Linux (Fedora/CentOS/RHEL):**

    ```bash
    sudo dnf install sshfs
    ```

*   **macOS:**

    1.  Install Homebrew (if you don't have it):  [https://brew.sh/](https://brew.sh/)
    2.  Install osxfuse and sshfs:

        ```bash
        brew install osxfuse
        brew install sshfs
        ```
        *Note:*  After installing osxfuse on macOS, you may need to go to System Preferences -> Security & Privacy and allow the osxfuse software to load.

*   **Windows:**

    *   SSHFS is more complex to set up on Windows. Consider using a graphical SFTP client like FileZilla or WinSCP for file transfer.  If you still prefer SSHFS, you can try using a pre-built SSHFS implementation for Windows, or consider using the Windows Subsystem for Linux (WSL) and installing SSHFS within your Linux distribution.

**3. Establishing the SSHFS Connection:**

Once SSHFS is installed, you can mount the HPCC filesystem. Here's the general command syntax:

```bash
sshfs <username>@<hpcc_head_node>:<remote_directory> <local_mount_point> [options]
```

*   `<username>`:  Your HPCC username.
*   `<hpcc_head_node>`: The address of the HPCC head node. (You can typically access this through the web console or by using `ursa-major.ucr.edu`)
*   `<remote_directory>`: The path to the directory on the HPCC you want to mount (e.g., `/home/<username>`).
*   `<local_mount_point>`: The path to the local directory on your laptop where you want to mount the filesystem (e.g., `~/hpcc_mount`).
*   `[options]` (Optional): Additional SSHFS options (explained below).

**Example:**

Let's say your HPCC username is `johndoe`, you want to connect to the HPCC head node (`ursa-major.ucr.edu`), you want to mount your home directory (`/home/johndoe`), and your local mount point is `~/hpcc_mount`.  The command would be:

```bash
sshfs johndoe@ursa-major.ucr.edu:/home/johndoe ~/hpcc_mount
```

**Common SSHFS Options:**

*   `-o allow_other`:  Allows other users on your system to access the mounted filesystem.  Use with caution, as this can have security implications.
*   `-o reconnect`:  Attempts to reconnect if the connection is lost.
*   `-o idmap=user`:  Maps the remote user ID to your local user ID.  This can help avoid permission issues.
*   `-o cache=yes`: Enables client-side caching.

**Example with Options:**

```bash
sshfs johndoe@ursa-major.ucr.edu:/home/johndoe ~/hpcc_mount -o allow_other,reconnect,idmap=user
```

**Security Note Regarding Passwords:**

The above commands will prompt you for your HPCC password. To avoid this, you can set up SSH key-based authentication. This is the recommended approach for security reasons.  Refer to the HPCC documentation or research-computing@ucr.edu to learn how to set up SSH key-based authentication.  Once set up, SSHFS will automatically use your SSH key to authenticate.

**4. Using the Mounted Filesystem:**

Once the connection is established, you can access the files and directories on the HPCC through your local mount point (`~/hpcc_mount` in the example). You can use your normal file manager or command-line tools to navigate, create, edit, and delete files.

**Example:**

```bash
cd ~/hpcc_mount
ls -l
mkdir new_directory
cp local_file.txt remote_file.txt
```

**5. Disconnecting the SSHFS Connection:**

When you're finished, you should unmount the filesystem. This is important to prevent data corruption and free up resources.

```bash
umount ~/hpcc_mount
```

**Important:** Make sure you are *not* in the mounted directory when you run the `umount` command.  Change to another directory first (e.g., `cd ~`).

**6. Troubleshooting:**

*   **Permission Denied:**  Ensure you have the correct permissions on the HPCC filesystem.  Use `chmod` on the HPCC if necessary. The `-o idmap=user` option can sometimes resolve permission issues.
*   **Connection Refused:** Verify that the HPCC head node address is correct and that the cluster is accessible. Check your network connection.
*   **"Device is busy" Error:**  Make sure no files or processes are accessing the mounted filesystem. Close any open files and change out of the mounted directory before unmounting.
*   **Slow Performance:** SSHFS can be slower than directly accessing the HPCC.  The `-o cache=yes` option can help improve performance. Consider using `rsync` or `scp` for large file transfers.
*   **Mount point doesn't exist:** Ensure that the local directory you are trying to mount to already exists

**7. Conclusion & Additional Resources:**

SSHFS provides a convenient way to access and manage files on the UCR HPCC from your local laptop.  By following these instructions, you can establish a secure connection and work with your data seamlessly.  Remember to disconnect the connection when you're finished.

**Additional Resources:**

*   **UCR Research Computing Documentation:** (Contact research-computing@ucr.edu) for documentation specific to the UCR HPCC.
*   **SSHFS Man Page:** `man sshfs` (for detailed information about SSHFS options).
*   **HPCC Documentation**: (Contact research-computing@ucr.edu) for documentation and guides on accessing and using the HPCC cluster.
*   **UCR Research Computing Slack:** [https://ucr-research-compute.slack.com/](https://ucr-research-compute.slack.com/)

For further assistance or questions, please contact UCR Research Computing at research-computing@ucr.edu.
