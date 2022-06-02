using System.Collections;
using System.Collections.Generic;
using UnityEngine;

using System;
using System.Net;
using System.Net.Sockets;
using System.Threading;

public class UDPListener : MonoBehaviour
{
    Thread receiveThread;
    UdpClient client;
    public int port;

    private Queue m_MsgQueue;

    // Start is called before the first frame update
    void Start()
    {
        Debug.Log("Starting UDP Listener");
        m_MsgQueue = new Queue();
        receiveThread = new Thread(new ThreadStart(SocketReceive));
        receiveThread.IsBackground = true;
        receiveThread.Start();
    }

    // Update is called once per frame
    void Update()
    {
        if (m_MsgQueue.Count != 0){
            byte[] myData = (byte[])m_MsgQueue.Dequeue();
            float[] myFloats = DecodeData(myData);
            // Convert to right hand coordinates
            transform.position = new Vector3(myFloats[0], myFloats[1], -myFloats[2]);
            float[] eul = {myFloats[3], myFloats[4], myFloats[5]};
            transform.rotation = ZXZEulerToQuaternion(eul); // Convert z - x' - z'' euler angles to a quaternion
        }   
    }

    private void SocketReceive()
    {
        client = new UdpClient(port);
        while (true)
        {
            try
            {
                IPEndPoint anyIP = new IPEndPoint(IPAddress.Any, 0);
                byte[] m_RecvData = client.Receive(ref anyIP);
                m_MsgQueue.Enqueue(m_RecvData);
            }
            catch (Exception err)
            {
                Console.WriteLine(err.ToString());
            }
        }
    }

    private float[] DecodeData(byte[] data)
    {
        if (data.Length != 4*6) // Check if data is 6 of 4-byte floats
        {
            Debug.Log("Bad data");
            return new float[6];
        }
        float[] returnData = new float[data.Length / 4];
        Buffer.BlockCopy(data, 0, returnData, 0, data.Length);
        for (int i = 0; i<3; i++)
        {
            returnData[i] = returnData[i];
        }
        return returnData;
    }


    // Adapted from NASA-TM-74839
    // Y and Z swapped due to left-hand coordinate system in unity
    private Quaternion ZXZEulerToQuaternion(float[] eul)
    {
        eul[1] = eul[1] + (float)Math.PI/2.0f;
        float w = (float)(Math.Cos(0.5f*eul[2])*Math.Cos(0.5f*(eul[0] + eul[1])));
        float x = (float)(Math.Sin(0.5f*eul[2])*Math.Cos(0.5f*(eul[0] - eul[1])));
        float y = (float)(Math.Sin(0.5f*eul[2])*Math.Sin(0.5f*(eul[0] - eul[1])));
        float z = (float)(Math.Cos(0.5f*eul[2])*Math.Sin(0.5f*(eul[0] + eul[1])));
        return new Quaternion(x, y, z, w);
    }
}