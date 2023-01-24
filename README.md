# soundloc_kmeans
Lightweight multiple sound source localization, based on a triangular microphone array.

Based on a combination of the [original soundloc package](https://github.com/balkce/soundloc) and Luis Gato's [kmeans-based improvement](https://github.com/lmiguelgato/DAP_project).

Uses JACK for input/audio audio server.

Outputs localizations through the JsonMsg topic, described by the json_msgs node.

Configured by soundloc_config.yaml:
* distance between microphones
* maximum number of sources
* graphical representation of localizations
* energy threshold to consider change of inactive and active states
* energy threshold to trigger activity
* coherence threshold (in degrees) between local microphone-pair estimations
* automatic connection to JACK inputs


## Dependencies
Packages that can be installed trough apt official repositories:
* libjack-jackd2-dev: JACK development libraries
* libfftw3-dev: a very fast FFT C/C++ implementation
* libjsoncpp-dev: for constructing JSON structures

