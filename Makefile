coverage: 
	/Applications/Julia-1.7.app/Contents/Resources/julia/bin/julia -e 'using Pkg; Pkg.activate("."); Pkg.test(;test_args=[""],coverage=true)'
	/Applications/Julia-1.7.app/Contents/Resources/julia/bin/julia -e 'using Coverage; coverage=process_folder(); LCOV.writefile("lcov.info", coverage); clean_folder(".")'
	