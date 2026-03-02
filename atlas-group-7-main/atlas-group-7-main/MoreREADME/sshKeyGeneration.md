In this file we will create an ssh-key to be used for with [git](https://git-scm.com/video/what-is-git)

1.  Open a terminal in your VM (or if you using your system/FP2 computer, terminal inthose)
2.  Run the command (you can replace your field `"your_email@example.com"` with your email or something that will contain your name ):
    ```
    ssh-keygen -t rsa -b 4096 -C "your_email@example.com" -f ~/.ssh/id_rsa.github
    ```
    **Keep the `-f ~/.ssh/id_rsa.github` as we need this particular name for the utput file for your git to work without any additional fixes** 
    (for people interested see `~/.ssh/config` file).

    * Do not use any password (type enter twice)
    * You should see something like this:
        ```
        Generating public/private rsa key pair.
        Enter passphrase (empty for no passphrase): 
        Enter same passphrase again: 
        Your identification has been saved in /home/fp2student/.ssh/id_rsa.github.
        Your public key has been saved in /home/fp2student/.ssh/id_rsa.github.pub.
        The key fingerprint is:
        SHA256:vA1dpk4MRPhUgN177EtE7NksMi7MYP4/FzAmeDTjsf0 IamFP2student@example.com
        The key's randomart image is:
        +---[RSA 4096]----+
        |       =++..     |
        |      o.B . o    |
        |       *.* =o+   |
        |      +.*+O+B o  |
        |     o =S+=X .   |
        |      . +*. E    |
        |       ...o. o   |
        |        . . o    |
        |         ..o     |
        +----[SHA256]-----+
        ```
    * The key you will get will be different from this example so the actual 'image' you get will be different
3. Type `cat ~/.ssh/id_rsa.github.pub` to see your public key. You should see something like this:
    ```
    ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAACAQDTBH3y33eZC5nwMKiji57x0GHS7LYYPgZ7LcYeN7Va+P37vXHmmnC\
    Kva1vFiibOIZIInyQGSk28876ks5F9osSInNUPTqzQKVhK+yLD+iNNYmvSFCnG75C6+knCE1ibpPI/4WsYkptjfd4+m\
    cJIMu7EaXAP5YUn9VMMhSxsDMy3ZrwxCd9jUwUTJUbPdHszodoOJHMWZLGFnj8p5CWxQv2ViAt/6SEn+V1u5sKVkJsJ\
    F9dP6urpsbWUrtub596j+7UCio2EqXtFgsi8SrpSjqnocjVK7xE+2dtj2z/P8yl5T+RInD8WZVmC0hsGcb+fM5hOA1p\
    py0p6WO2jqyco1X6pquTWdrQtGoQg+lXi/0JF3PuEoHJPw07km3j7NOJVZxCQ2kStOI8h3RtvCQzq/RuMjsTzXmmhMP\
    NL1MiFGQMz1mwsbAVnqX+O50mLVfTj7AUj6+8JcJYquUbcQmxb7X8VUALTy+/eFH0zzccTbyXYAAVvrujorCWWomOav\
    wQuIaaGFH9q75frBpYcJ7KBax1/KtI4rlc2ZMWyEenE/xXnsGbeGHXxew7DYYUFOx/O429QvAOOcJDnnnopcxNAqblq\
    QHNzVm/eZ8VK6GQ7kDiqQyPlPMSYQiK4Ul6+E87VfBuVT0szWdIKqoKAk/c7IUzNOX8ZK3juFm1NupCndRbBQ== Iam\
    FP2student@example.com
    ```
3. Open your github account
4. Get to settings of your github account (upper-rght corner, next to your account icon)
5. Go to `SSH and GPG keys` tab
6. Click `New SSH key` button
    * Choose any name you want for your key
    * In the field `Key` paste the output you got from `cat ~/.ssh/id_rsa.github.pub` command in point 3. 
    * Click `Add SSH key` to finish the procedure
7. Now you should have the ssh-key setup to run with your git! Ask tutor for ways to verify that this is indeed the case!