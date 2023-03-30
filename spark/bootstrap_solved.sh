#!/bin/bash
set -x
 
cat > /var/tmp/fix-bootstap.sh <<'EOF'
#!/bin/bash
set -x
 
while true; do
    NODEPROVISIONSTATE=`sed -n '/localInstance [{]/,/[}]/{
    /nodeProvisionCheckinRecord [{]/,/[}]/ {
    /status: / { p }
    /[}]/a
    }
    /[}]/a
    }' /emr/instance-controller/lib/info/job-flow-state.txt | awk ' { print $2 }'`
 
    if [ "$NODEPROVISIONSTATE" == "SUCCESSFUL" ]; then
        echo "Running my post provision bootstrap"
        # Enter your code here
        sudo yum install python3-devel   # for python3.x installs
        sudo yum install python3.8-devel   # for python3.x installs
        sudo python3 -m pip install --upgrade pip
        sudo python3 -m pip install sklearn
        sudo python3 -m pip install numpy
        sudo python3 -m pip install pandas
        sudo python3 -m pip install scikit-build
        sudo python3 -m pip install geopandas
        sudo python3 -m pip install pyrosm
        sudo python3 -m pip install pandas geopy pyarrow shapely geopandas rtree matplotlib mapclassify scipy fsspec s3fs openpyxl boto3 statsmodels h3 h3-pyspark h3pandas osmium osmnx
        echo '-------BOOTSTRAP COMPLETE-------' 
 
        exit
    else
        echo "Sleeping Till Node is Provisioned"
        sleep 10
    fi
done
 
EOF
 
chmod +x /var/tmp/fix-bootstap.sh
nohup /var/tmp/fix-bootstap.sh  2>&1 &

