This code is an update version of https://github.com/yujieqing/SegmentsComputation.

The code is for geometric (or segment) matrix computation in voxel-based ionosphere tomography inversion. Both VIST and AP method were implemented in this code. 

**AP approach**: see, Hong J, Kim Y H, Chung J K, et al. (2017) Tomography reconstruction of ionospheric electron density with empirical orthonormal functions using Korea GNSS network. J Astron Space Sci 34(1): 7-17.

**VIST approach**: introduced in the manuscirpt named "Fast determination of geometric matrix in ionosphere tomographic inversion with unevenly spaced curvilinear voxels" submitted to GPS Solutions.

Please use VS2019 to compile it.

The output can be visualized by paraview or VisIT software after conversion using scripts 'toVtp.py'. To use the script, python3.7 is required and vtk libary should be installed.  



	/main.cpp: the main function
	/Common.h,Common.cpp: common functions for both methods
	/AP.h,tradition.cpp: AP method. 
	/VIST.h,VIST.cpp: VIST method
	/toVtp.py, vtkXMLFile.py: scripts for converting the resultant segments into .vtp format, which can be recognized by paraview(https://www.paraview.org/) or VisIT (https://wci.llnl.gov/simulation/computer-codes/visit/) software.
	/input
		/Conf.txt：file that includes input parameters
		/rays: files that contain the path coordinates between satellite and ground station.
	/output
		time_cost.txt: time cost for both two methods
		AP_sparse(i).txt: the output of AP method in sparse matrix format, where ray ID, voxel ID, and the segment length are recorded.
		VIST_sparse(i).txt: the output of VIST method in sparse matrix format
		AP_centSize(i).txt: the output of AP method in centSize format, which can be read and converted into .vtp by toVtp.py
		VIST_centSize(i).txt: the output of VIST method in centSize format, which can be read and converted into .vtp by toVtp.py
		
For any question, please contact Dr. Jieqing YU (yujieqing@cumt.edu.cn)
