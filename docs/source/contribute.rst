Contribute to OpenMSCG
======================

To be able to contribute to the OpenMSCG community by submitting new features or
fixing bugs, you need to enable the GitLab account access described in this page.

Setup SSH Key on GitLab
-----------------------

Login into the GitLab portal (`Link <https://software.rcc.uchicago.edu/git/>`_)
from the RCC Software website. If you don't have a GitLab account, please email 
`Yuxing <yuxing@uchicago.edu>`_ to get one for free.

Setup your SSH public key to connect to GitLab (`instructions <https://help.github.
com/en/articles/adding-a-new-ssh-key-to-your-github-account>`_). If you are unfamiliar
with SSH key pairs or how to create them, please read the instructions
`here <https://www.digitalocean.com/docs/droplets/how-to/add-ssh-keys/create-
with-openssh/>`_.

Check that the SSH key was paired successfully::

    ssh -T git@software.rcc.uchicago.edu
    Welcome to GitLab, @yuxing!

If the SSH key is working, a welcome message from GitLab will be shown.
