# Wind-Aware AirSim

**WARNING** These modifications are a work in progress, and should be treated as such. They may have bugs and have not been extensively tested. The entire setup has been tested using Ubuntu 18.04 and 20.04 and Unreal Engine 4.25. The most up to date version of AirSim uses Unreal Engine 4.27 - we plan to update and rebase off of the latest changes to AirSim in the future. Due to the WIP nature and relative instability of this project, no compiled binaries are available. Developer access to Unreal Engine's source (easily available) is necessary.

## 1 Installing Unreal Engine

1. Clone Unreal Engine (version 4.25)

```shell
git clone -b 4.25 git@github.com:EpicGames/UnrealEngine.git
```

2. Build Unreal (this will take awhile)

```shell
cd UnrealEngine
./Setup.sh
./GenerateProjectFiles.sh
make
```

These steps are identical to those provided by the AirSim authors themselves, visible [here](https://microsoft.github.io/AirSim/build_linux/).

## 2 Installing Wind-Aware AirSim

1. Clone Wind-Aware AirSim

```shell
git clone -b spatial-temporal-wind git@github.com:CoRAL-OSU/AirSim_NRI.git
```

2. Build Wind-Aware AirSim

```shell
cd AirSim_NRI
./setup.sh
./build.sh
```

## 3 Installing Wind-Aware PX4 Flight Controller


1. Clone Wind-Aware PX4

```shell
git clone -b wind_aware git@github.com:CoRAL-OSU/wind_aware_px4.git --recursive
```

2. Build Wind-Aware PX4

```shell
cd wind_aware_px4
make px4_sitl_default none_iris
```

Further information on PX4 SITL is available from AirSim's docs, [here](https://microsoft.github.io/AirSim/px4_sitl/).

## 4 Installing Wind-Aware QGroundControl

1. Clone Wind-Aware QGC

```shell
git clone -b wind-display-ui git@github.com:CoRAL-OSU/wind_aware_QGC.git --recursive
```

2. Build Wind-Aware QGC

This must be done using Qt Creator right now. This can be done using a process identical to what is provided by QGC's authors, available [here](https://dev.qgroundcontrol.com/master/en/getting_started/index.html). Our version of Wind-Aware QGC is for Qt version 5.15.2


# Using Wind-Aware AirSim

## Basics

1. Launch the UE4 editor found at the path `<unreal install dir>/Engine/Binaries/Linux/UE4Editor`
2. Browse for the "Blocks" environment included with AirSim, found under `<airsim install dir>/Unreal/Environments/Blocks`
3. Convert/build the environment if necessary. The "convert-in-place" option may be under a "more options" tab.
4. Once inside the editor, an example map of Oklahoma State University is provided. This map is generated using data from OpenStreetMap, and is not entirely collision-friendly. Information on building your own environment is found [here](https://microsoft.github.io/AirSim/unreal_custenv/). Alternatively, the "FlyingExampleMap.umap" level can be selected, for a small test environment.
5. Press the "Play" option to begin the simulation. This uses settings found in AirSim's `settings.json` file.

## Wind Configuration

Wind settings are configured in the `settings.json` file. An example settings file is provided in the root of the project directory. An example configuration is shown below:

```
"Wind": {
  "EnableLog": true,
  "WindDataPath": "/home/coral/winddata/wind_001.txt",
  "UpdatePeriod": 2.5,
  "DefaultWind": {
    "X": 5,
    "Y": 4,
    "Z": 3
  }
}
```

`WindDataPath` should point to where the spatio-temporally varying wind data is located. `UpdatePeriod` is how many seconds the simulator waits before loading the next data file. `DefaultWind` is used if the aircraft is outside the provided area, or if no wind data is provided.

The simulator expects temporal data to advance with a specific naming scheme. Subsequent data files should be numbered in the manner "wind_001.txt", "wind_002.txt" ... "wind_013.txt", etc. All data should use the same XYZ bounds and step sizes. The wind data format is adapted [here](https://github.com/ethz-asl/rotors_simulator/wiki/Adding-a-custom-wind-field-to-your-world). Example wind file may be requested by email @asma.tabassum@okstate.edu.

If you are using this simulator within the research for your publication, please cite:
```
@inproceedings{tabassum2022preliminary,
  title={Preliminary Design of Wind-Aware sUAS Simulation Pipeline for Urban Air Mobility},
  author={Tabassum, Asma and DeSantis, Max and Bai, He and Fala, Nicoletta},
  booktitle={AIAA AVIATION 2022 Forum},
  pages={3872},
  year={2022}
}
```

---

# Welcome to AirSim

AirSim is a simulator for drones, cars and more, built on [Unreal Engine](https://www.unrealengine.com/) (we now also have an experimental [Unity](https://unity3d.com/) release). It is open-source, cross platform, and supports software-in-the-loop simulation with popular flight controllers such as PX4 & ArduPilot and hardware-in-loop with PX4 for physically and visually realistic simulations. It is developed as an Unreal plugin that can simply be dropped into any Unreal environment. Similarly, we have an experimental release for a Unity plugin.

Our goal is to develop AirSim as a platform for AI research to experiment with deep learning, computer vision and reinforcement learning algorithms for autonomous vehicles. For this purpose, AirSim also exposes APIs to retrieve data and control vehicles in a platform independent way.

**Check out the quick 1.5 minute demo**

Drones in AirSim

[![AirSim Drone Demo Video](docs/images/demo_video.png)](https://youtu.be/-WfTr1-OBGQ)

Cars in AirSim

[![AirSim Car Demo Video](docs/images/car_demo_video.png)](https://youtu.be/gnz1X3UNM5Y)


## How to Get It

### Windows
[![Build Status](https://github.com/microsoft/AirSim/actions/workflows/test_windows.yml/badge.svg)](https://github.com/microsoft/AirSim/actions/workflows/test_windows.yml)
* [Download binaries](https://github.com/Microsoft/AirSim/releases)
* [Build it](https://microsoft.github.io/AirSim/build_windows)

### Linux
[![Build Status](https://github.com/microsoft/AirSim/actions/workflows/test_ubuntu.yml/badge.svg)](https://github.com/microsoft/AirSim/actions/workflows/test_ubuntu.yml)
* [Download binaries](https://github.com/Microsoft/AirSim/releases)
* [Build it](https://microsoft.github.io/AirSim/build_linux)

### macOS
[![Build Status](https://github.com/microsoft/AirSim/actions/workflows/test_macos.yml/badge.svg)](https://github.com/microsoft/AirSim/actions/workflows/test_macos.yml)
* [Build it](https://microsoft.github.io/AirSim/build_linux)

For more details, see the [use precompiled binaries](docs/use_precompiled.md) document. 

## How to Use It

### Documentation

View our [detailed documentation](https://microsoft.github.io/AirSim/) on all aspects of AirSim.

### Manual drive

If you have remote control (RC) as shown below, you can manually control the drone in the simulator. For cars, you can use arrow keys to drive manually.

[More details](https://microsoft.github.io/AirSim/remote_control)

![record screenshot](docs/images/AirSimDroneManual.gif)

![record screenshot](docs/images/AirSimCarManual.gif)


### Programmatic control

AirSim exposes APIs so you can interact with the vehicle in the simulation programmatically. You can use these APIs to retrieve images, get state, control the vehicle and so on. The APIs are exposed through the RPC, and are accessible via a variety of languages, including C++, Python, C# and Java.

These APIs are also available as part of a separate, independent cross-platform library, so you can deploy them on a companion computer on your vehicle. This way you can write and test your code in the simulator, and later execute it on the real vehicles. Transfer learning and related research is one of our focus areas.

Note that you can use [SimMode setting](https://microsoft.github.io/AirSim/settings#simmode) to specify the default vehicle or the new [ComputerVision mode](https://microsoft.github.io/AirSim/image_apis#computer-vision-mode-1) so you don't get prompted each time you start AirSim.

[More details](https://microsoft.github.io/AirSim/apis)

### Gathering training data

There are two ways you can generate training data from AirSim for deep learning. The easiest way is to simply press the record button in the lower right corner. This will start writing pose and images for each frame. The data logging code is pretty simple and you can modify it to your heart's content.

![record screenshot](docs/images/record_data.png)

A better way to generate training data exactly the way you want is by accessing the APIs. This allows you to be in full control of how, what, where and when you want to log data.

### Computer Vision mode

Yet another way to use AirSim is the so-called "Computer Vision" mode. In this mode, you don't have vehicles or physics. You can use the keyboard to move around the scene, or use APIs to position available cameras in any arbitrary pose, and collect images such as depth, disparity, surface normals or object segmentation.

[More details](https://microsoft.github.io/AirSim/image_apis)

### Weather Effects

Press F10 to see various options available for weather effects. You can also control the weather using [APIs](https://microsoft.github.io/AirSim/apis#weather-apis). Press F1 to see other options available.

![record screenshot](docs/images/weather_menu.png)

## Tutorials

- [Video - Setting up AirSim with Pixhawk Tutorial](https://youtu.be/1oY8Qu5maQQ) by Chris Lovett
- [Video - Using AirSim with Pixhawk Tutorial](https://youtu.be/HNWdYrtw3f0) by Chris Lovett
- [Video - Using off-the-self environments with AirSim](https://www.youtube.com/watch?v=y09VbdQWvQY) by Jim Piavis
- [Webinar - Harnessing high-fidelity simulation for autonomous systems](https://note.microsoft.com/MSR-Webinar-AirSim-Registration-On-Demand.html) by Sai Vemprala
- [Reinforcement Learning with AirSim](https://microsoft.github.io/AirSim/reinforcement_learning) by Ashish Kapoor
- [The Autonomous Driving Cookbook](https://aka.ms/AutonomousDrivingCookbook) by Microsoft Deep Learning and Robotics Garage Chapter
- [Using TensorFlow for simple collision avoidance](https://github.com/simondlevy/AirSimTensorFlow) by Simon Levy and WLU team

## Participate

### Paper

More technical details are available in [AirSim paper (FSR 2017 Conference)](https://arxiv.org/abs/1705.05065). Please cite this as:
```
@inproceedings{airsim2017fsr,
  author = {Shital Shah and Debadeepta Dey and Chris Lovett and Ashish Kapoor},
  title = {AirSim: High-Fidelity Visual and Physical Simulation for Autonomous Vehicles},
  year = {2017},
  booktitle = {Field and Service Robotics},
  eprint = {arXiv:1705.05065},
  url = {https://arxiv.org/abs/1705.05065}
}
```

### Contribute

Please take a look at [open issues](https://github.com/microsoft/airsim/issues) if you are looking for areas to contribute to.

* [More on AirSim design](https://microsoft.github.io/AirSim/design)
* [More on code structure](https://microsoft.github.io/AirSim/code_structure)
* [Contribution Guidelines](CONTRIBUTING.md)

### Who is Using AirSim?

We are maintaining a [list](https://microsoft.github.io/AirSim/who_is_using) of a few projects, people and groups that we are aware of. If you would like to be featured in this list please [make a request here](https://github.com/microsoft/airsim/issues).

## Contact

Join our [GitHub Discussions group](https://github.com/microsoft/AirSim/discussions) to stay up to date or ask any questions.

We also have an AirSim group on [Facebook](https://www.facebook.com/groups/1225832467530667/). 


## What's New

- [Python wrapper for Open AI gym interfaces.](https://github.com/microsoft/AirSim/pull/3215)
- [Python wrapper for Event camera simulation](https://github.com/microsoft/AirSim/pull/3202)
- [Voxel grid construction](https://github.com/microsoft/AirSim/pull/3209)
- [Programmable camera distortion](https://github.com/microsoft/AirSim/pull/3039)
- [Wind simulation](https://github.com/microsoft/AirSim/pull/2867)
- [Azure development environment with documentation](https://github.com/microsoft/AirSim/pull/2816)
- ROS wrapper for [multirotor](https://github.com/microsoft/AirSim/blob/master/docs/airsim_ros_pkgs.md) and [car](https://github.com/microsoft/AirSim/pull/2743).

For complete list of changes, view our [Changelog](docs/CHANGELOG.md)

## FAQ

If you run into problems, check the [FAQ](https://microsoft.github.io/AirSim/faq) and feel free to post issues in the  [AirSim](https://github.com/Microsoft/AirSim/issues) repository.

## Code of Conduct

This project has adopted the [Microsoft Open Source Code of Conduct](https://opensource.microsoft.com/codeofconduct/). For more information see the [Code of Conduct FAQ](https://opensource.microsoft.com/codeofconduct/faq/) or contact [opencode@microsoft.com](mailto:opencode@microsoft.com) with any additional questions or comments.


## License

This project is released under the MIT License. Please review the [License file](LICENSE) for more details.


