import socket
import struct
import multiprocessing
import datetime

# Listen for data from Dirt/F1 over specified UDP port
# Extract float values, write to CSV
# Data format: https://docs.google.com/spreadsheets/d/1eA518KHFowYw7tSMa-NxIFYpiWe5JXgVVQ_IMs7BVW0/edit#gid=0

UDP_PORT = 20777

def log(client):
    dt = datetime.datetime.now()
    datetimestr = dt.strftime("%Y%m%dT%H%M%S")
    
    with open("logs/" + datetimestr + ".log", "x") as f:
        while True:
            try:
                packet = client.recv(1024)
                data = struct.unpack_from("<63f", packet, 0)
                csvline = ""
                for line in data:
                    csvline += str(line) + ","
                csvline, _ = csvline.rsplit(",", maxsplit=1) # remove last comma
                f.write(csvline + "\n") # Write to file
                f.flush()
            except Exception as e:
                pass # Ignore errors

if __name__ == "__main__":
    addr = ("", UDP_PORT)
    client = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    client.bind(addr)

    p1 = multiprocessing.Process(target=log, args=(client, ))
    p1.start()

    input("Press enter to stop") # Wait for logging stop
    p1.kill()
    client.close()
