if [ "$#" -lt 2 ]; then
    echo "Usage "$0" <databasename> <h5files>"
else
echo "Last update in database "$(cat $1 | sort -k 5 | tail -n 1 | awk -F' ' '{ print $5" "$6 }')
echo "Timestamps in h5 files:"
. h5toObsstart.sh ${@:2}
echo "If time stamp in files are newer than the database, please download a new version:
wget --no-check-certificate https://proxy.lofar.eu/array_status/STATIONS/DATA/hardware_states_latest.txt -O DATA/hardware_states_latest.txt
If this doesn't fix it, request a new database by sending in a support ticket to SOS."
fi
