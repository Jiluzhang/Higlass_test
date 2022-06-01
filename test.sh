git -C /usr/local/Homebrew/Library/Taps/homebrew/homebrew-core fetch --unshallow
#remote: Enumerating objects: 1170676, done.
#remote: Counting objects: 100% (1170246/1170246), done.
#remote: Compressing objects: 100% (362598/362598), done.
#remote: Total 1161925 (delta 806276), reused 1151960 (delta 796382), pack-reused 0
#Receiving objects: 100% (1161925/1161925), 451.97 MiB | 10.12 MiB/s, done.
#Resolving deltas: 100% (806276/806276), completed with 6817 local objects.
#From https://github.com/Homebrew/homebrew-core
#   683982a0120..f95e28707f9  master     -> origin/master

# find ip: https://ipaddress.com/website/github.com
sudo vi /private/etc/hosts
# 140.82.112.4  github.com

brew update
#==> Downloading https://ghcr.io/v2/homebrew/portable-ruby/portable-ruby/blobs/sha256:0cb1cc7af109437fe0e020c9f3b7b95c3c709b140bde9f991ad2c1433496dd42
######################################################################### 100.0%
#==> Pouring portable-ruby-2.6.8.yosemite.bottle.tar.gz
#Updated 1 tap (homebrew/core).
#==> Updated Formulae
#gopass
#
#You have 18 outdated formulae installed.
#You can upgrade them with brew upgrade
#or list them with brew outdated.

brew install pkg-config

pip install --upgrade setuptools

pip install higlass-manage  # finally ok!!!! ^-^
