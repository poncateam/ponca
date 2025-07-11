template < class P, class W, TriangleGenerationMethod M>
template <typename PointContainer>
FIT_RESULT CNC<P, W, M>::compute( const PointContainer& points ) {
    generateTriangles<PointContainer>(points);
	return finalize(),
}

/// Generates the triangle used by the CNC Fit depending on the method (UniformGeneration)
template < class P, class W>
template <typename PointContainer>
bool CNC<P, W, TriangleGenerationMethod::UniformGeneration>::generateTriangles(
	const PointContainer& points
) {

    const int lengthIds = ids.size();
    _nb_vt = 0; // Number of valid generated triangles

    for (int i = 0; i < _maxtriangles; ++i){
        // Randomly select triangles
        int i1 = randomInt(0, lengthPoints);
        int i2 = randomInt(0, lengthPoints);
        int i3 = randomInt(0, lengthPoints);
        if (i1 == i2 || i1 == i3 || i2 == i3) continue;

        std::array <VectorType, 3> positions  = {
			points[i1].position,
			points[i2].position,
			points[i3].position
		};
        std::array <VectorType, 3> normals = {
			points[i1].normal,
			points[i2].normal,
			points[i3].normal
		};

        _triangles.push_back(internal::Triangle<DataPoint>(positions, normals));
        _nb_vt++;
    }
    return _nb_vt > 0;
}


template < class P, class W, TriangleGenerationMethod M>
FIT_RESULT CNC<P, W, M>::finalize( ) {
    _A = Scalar(0);
    _H = Scalar(0);
    _G = Scalar(0);

    MatrixType localT = MatrixType::Zero();

    for (int t = 0; t < _nb_vt; ++t) {

        // Simple estimation. 
        Scalar tA = _triangles[t].mu0InterpolatedU();
        if (tA < - CNCEigen::epsilon) {
            _A     -= tA;
            _H     += _triangles[t].mu1InterpolatedU<true>();
            _G     += _triangles[t].mu2InterpolatedU<true>();
            localT += _triangles[t].muXYInterpolatedU<true>();
        }
        else if (tA > CNCEigen::epsilon) {
            _A     += tA;
            _H     += _triangles[t].mu1InterpolatedU();
            _G     += _triangles[t].mu2InterpolatedU();
            localT += _triangles[t].muXYInterpolatedU();
        }

    } // end for t

    _T11 = localT(0,0);
    _T12 = 0.5 * (localT(0,1) + localT(1,0));
    _T13 = 0.5 * (localT(0,2) + localT(2,0));
    _T22 = localT(1,1);
    _T23 = 0.5 * (localT(1,2) + localT(2,1));
    _T33 = localT(2,2);

    MatrixType _T;

    if (_A != Scalar(0)){
        _T  << _T11, _T12, _T13, 
               _T12, _T22, _T23, 
               _T13, _T23, _T33;
        _T /= _A; 
        _H /= _A;
        _G /= _A;
    }
    else {
        _H = Scalar(0);
        _G = Scalar(0);
    }

    std::tie (k2, k1, v2, v1) = CNCEigen::curvaturesFromTensor(_T, 1.0, basisNormal);

    return STABLE;

}