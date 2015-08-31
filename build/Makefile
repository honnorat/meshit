cmake-release:
	(mkdir -p Release ; cd Release ; cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_CXX_FLAGS="-O3" ../..)

cmake-debug:
	(mkdir -p Debug ; cd Debug ; cmake -DCMAKE_BUILD_TYPE=DEBUG ../..)

cmake-profile:
	(mkdir -p Profile ; cd Profile ; cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_CXX_FLAGS="-O3 -g" ../..)

make-release: cmake-release
	(cd Release ; $(MAKE) -j8 -f Makefile)

make-profile: cmake-profile
	(cd Profile ; $(MAKE) -j8 -f Makefile)

make-debug: cmake-debug
	(cd Debug ; $(MAKE) -j8 -f Makefile)

clean-release:
	(cd Release ; $(MAKE) clean)

clean-debug:
	(cd Debug ; $(MAKE) clean)

clean: clean-debug
