a=`ps -ef | grep 19959 | grep -v grep`

if [ ! "$a" ]; then

    echo "Start Tunnel"
    ssh -fN -R 19959:localhost:22 master@airglow.ece.illinois.edu
fi

echo "Tunnel Open"
