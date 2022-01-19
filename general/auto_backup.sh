#!/bin/sh
# A personal script to run rsync over SSH for specific folders that should be 
# the same for both my personal laptop, and my VUW mac. 
#
# useful references:
# + anyconnect on cli
# https://superuser.com/questions/649614/connect-using-anyconnect-from-command-line
# + osx passwords on keychain using cli
# https://www.netmeister.org/blog/keychain-passwords.html
# + accessing passwords from keychain using cli
# https://macromates.com/blog/2006/keychain-access-from-shell/
#

# Specify login credentials
login="chowbr"
pword=$(security 2>&1 >/dev/null find-generic-password -ga ${login} | ruby -e 'print $1 if STDIN.gets =~ /^password: "(.*)"$/')

# VPN using anyconnect binary, redirect outputs to dev/null
echo "connecting to VPN"
printf "y\n${login}\n${pword}\ny" | /opt/cisco/anyconnect/bin/vpn -s connect vpn.vuw.ac.nz >/dev/null

# Rsync the necessary files based on hardcoded directory structure
echo "rsyncing"
LAP_SUB="/Users/Chow/Documents/academic/vuw"
VUW_SUB="/Users/chowbr/Documents/subduction"
SSH="chowbr@C02SX063GG7D.staff.vuw.ac.nz"

# rsync both ways so both computers are the same
for DIR in oficial
do
	# Laptop -> iMac
	# rsync -e ssh -avzp ${LAP_SUB}/${DIR} ${SSH}:${VUW_SUB}
	# iMac -> Laptop
	rsync -e ssh -avzp ${SSH}:${VUW_SUB}/${DIR} ${LAP_SUB}
done

echo "disconnecting from VPN"
/opt/cisco/anyconnect/bin/vpn disconnect >/dev/null
