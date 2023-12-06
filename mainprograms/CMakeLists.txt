add_executable(TestEigen mainEigenTest.cpp ${SOURCES} ${HEADERS})
target_link_libraries (TestEigen Eigen3::Eigen)

add_executable(Fem2021 main2021.cpp ${SOURCES} ${HEADERS})
target_link_libraries (Fem2021 Eigen3::Eigen)

add_executable(Integration mainintegration.cpp ${SOURCES} ${HEADERS})
target_link_libraries (Integration Eigen3::Eigen)

add_executable(Geometry maingeom.cpp ${SOURCES} ${HEADERS})
target_link_libraries (Geometry Eigen3::Eigen)

add_executable(Stiff stiff.cpp ${SOURCES} ${HEADERS})
target_link_libraries (Stiff Eigen3::Eigen)

add_executable(TestIntegrate TestIntegrate.cpp ${SOURCES} ${HEADERS})
target_link_libraries (TestIntegrate PRIVATE Eigen3::Eigen Catch2::Catch2)

add_executable(TestGeom TestGeom.cpp ${SOURCES} ${HEADERS})
target_link_libraries (TestGeom PRIVATE Eigen3::Eigen Catch2::Catch2)

add_executable(TestShape TestShape.cpp ${SOURCES} ${HEADERS})
target_link_libraries (TestShape PRIVATE Eigen3::Eigen Catch2::Catch2)

add_executable(TestOneDProblem TestOneDProblem.cpp ${SOURCES} ${HEADERS})
target_link_libraries (TestOneDProblem Eigen3::Eigen)

add_executable(TestTwoDProblem TestTwoDProblem.cpp ${SOURCES} ${HEADERS})
target_link_libraries (TestTwoDProblem Eigen3::Eigen)

add_executable(TestThreeDProblem TestThreeDProblem.cpp ${SOURCES} ${HEADERS})
target_link_libraries (TestThreeDProblem Eigen3::Eigen)

add_executable(ExerciseTwo ExerciseTwo.cpp ${SOURCES} ${HEADERS})
target_link_libraries (ExerciseTwo Eigen3::Eigen)

add_executable(ExerciseTwo ExerciseThree.cpp ${SOURCES} ${HEADERS})
target_link_libraries (ExerciseThree Eigen3::Eigen)
