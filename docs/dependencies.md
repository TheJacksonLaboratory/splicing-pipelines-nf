# Dependencies

## Rationale

Some pipelines require connection to external servers, to retrieve information, eg API requests to Ensembl.

In environments that have strict network permissions, eg access to the internet is blocked, this results in pipeline failure that is not related to the actual code but the network connection restrictions.

To check both explicit and non-explicit URL requests, [tcpdump](https://github.com/the-tcpdump-group/tcpdump) is a data-network packet analyzer computer program that runs under a command line interface. It allows the user to display TCP/IP and other packets being transmitted or received over a network to which the computer is attached.

tcpdump compiles and works on at least the following platforms:

- GNU/Linux
- {Mac} OS X / macOS
- Windows (requires WinPcap or Npcap, and Visual Studio with CMake)

### Install tcpdump

To install `tcpdump` run the following commands:

```bash
apt-get update && apt-get install -y tcpdump
```

> sudo might be required

### Launch tcpdump

To launch `tcpdump` before running nextflow open a new terminal and paste
the following command:

```bash
sudo tcpdump -i any -A -vv -s 0 | grep -Eo "> .*\.(ftp|https|http)" > connections.txt
```

All connections are being stored in `connections.txt`. This records the
URL for every packet transaction, so duplication are expected.

### Launch nextflow

In a new terminal launch your nextflow job. After it completes you can terminate the `tcpdump` process by pressing `ctrl + c`, and manually observe the connections file with:

```bash
cat connections.txt | grep -v net | grep -v ec2 | sort | uniq > unique_connections.txt
```

The resulting output should be (truncated):

```bash
> s3-3-w.amazonaws.com.https  <---
> s3-eu-west-1-r-w.amazonaws.com.https  <---
```

Highlighted with arrows (`<---`) are two connections to external servers, both `https` for amazon servers.
