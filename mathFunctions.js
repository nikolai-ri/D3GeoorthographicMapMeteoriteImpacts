// Credits to http://bl.ocks.org/ivyywang/7c94cb5a3accd9913263 for implementing these functions.
// I merely adapted to ES6 and created all the promises, to avoid problems with asynchronicity.

/***** ALL MATH FUNCTIONS ****/

const to_radians = Math.PI / 180;
const to_degrees = 180 / Math.PI;


// Helper function: cross product of two vectors v0&v1
function cross(v0, v1) {
    return [v0[1] * v1[2] - v0[2] * v1[1], 
    	    v0[2] * v1[0] - v0[0] * v1[2], 
    	    v0[0] * v1[1] - v0[1] * v1[0]];
}

//Helper function: dot product of two vectors v0&v1
function dot(v0, v1) {
	let sum = 0;
    for (let i = 0; v0.length > i; ++i) sum += v0[i] * v1[i];
    return sum;
}

// Helper function: 
// This function converts a [lon, lat] coordinates into a [x,y,z] coordinate 
// the [x, y, z] is Cartesian, with origin at lon/lat (0,0) center of the earth
function lonlat2xyz( coord ){
	return new Promise(function(res, rej){
		const lon = coord[0] * to_radians;
		const lat = coord[1] * to_radians;

		const x = Math.cos(lat) * Math.cos(lon);

		const y = Math.cos(lat) * Math.sin(lon);

		const z = Math.sin(lat);

		res([x, y, z]);
	});
}

// Helper function: 
// This function computes a quaternion representation for the rotation between to vectors
// https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Euler_angles_.E2.86.94_Quaternion
function quaternion(v0, v1) {
	return new Promise(function(res, rej){

		if(!v0 || !v1) rej;

	    let w = cross(v0, v1);  // vector pendicular to v0 & v1
	    const w_len = Math.sqrt(dot(w, w)); // length of w     

	    if (w_len == 0)
	    	return;

	    let theta = .5 * Math.acos(Math.max(-1, Math.min(1, dot(v0, v1)))),

	        qi  = w[2] * Math.sin(theta) / w_len; 
	        qj  = - w[1] * Math.sin(theta) / w_len; 
	        qk  = w[0]* Math.sin(theta) / w_len;
	        qr  = Math.cos(theta);

	    res(theta && [qr, qi, qj, qk]);
	});

}

// Helper function: 
// This functions converts euler angles to quaternion
// https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Euler_angles_.E2.86.94_Quaternion
function euler2quat(e) {
    return new Promise(function(res, rej){
    	
    	if(!e) rej;
    	let roll = .5 * e[0] * to_radians,
        pitch = .5 * e[1] * to_radians,
        yaw = .5 * e[2] * to_radians,

        sr = Math.sin(roll),
        cr = Math.cos(roll),
        sp = Math.sin(pitch),
        cp = Math.cos(pitch),
        sy = Math.sin(yaw),
        cy = Math.cos(yaw),

        qi = sr*cp*cy - cr*sp*sy,
        qj = cr*sp*cy + sr*cp*sy,
        qk = cr*cp*sy - sr*sp*cy,
        qr = cr*cp*cy + sr*sp*sy;

        res([qr, qi, qj, qk]);
    });
}

// This functions computes a quaternion multiply
// Geometrically, it means combining two quant rotations
// http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions/arithmetic/index.htm
function quatMultiply(q1, q2) {
	return new Promise(function(res, rej){
		if(!q1 || !q2) rej;

	    const a = q1[0],
	          b = q1[1],
	          c = q1[2],
	          d = q1[3],
	          e = q2[0],
	          f = q2[1],
	          g = q2[2],
	          h = q2[3];

	    res ([
	     a*e - b*f - c*g - d*h,
	     b*e + a*f + c*h - d*g,
	     a*g - b*h + c*e + d*f,
	     a*h + b*g - c*f + d*e]);

	});
}

// This function computes quaternion to euler angles
// https://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Euler_angles_.E2.86.94_Quaternion
function quat2euler(t){

	if(!t) return;

	return [ Math.atan2(2 * (t[0] * t[1] + t[2] * t[3]), 1 - 2 * (t[1] * t[1] + t[2] * t[2])) * to_degrees, 
			 Math.asin(Math.max(-1, Math.min(1, 2 * (t[0] * t[2] - t[3] * t[1])))) * to_degrees, 
			 Math.atan2(2 * (t[0] * t[3] + t[1] * t[2]), 1 - 2 * (t[2] * t[2] + t[3] * t[3])) * to_degrees
			]
}

/*  This function computes the euler angles when given two vectors, and a rotation
	This is really the only math function called with d3 code.

	v0 - starting pos in lon/lat, commonly obtained by projection.invert
	v1 - ending pos in lon/lat, commonly obtained by projection.invert
	o0 - the projection rotation in euler angles at starting pos (v0), commonly obtained by projection.rotate
*/

function eulerAngles(v0, v1, o0) {

	/*
		The math behind this:
		- first calculate the quaternion rotation between the two vectors, v0 & v1
		- then multiply this rotation onto the original rotation at v0
		- finally convert the resulted quat angle back to euler angles for d3 to rotate
	*/
	
	return new Promise(function (res, rej){
		let q1Res = 0,
		    q2Res = 0,
		    v0Res = 0,
		    v1Res = 0;
		
		euler2quat(o0).then(function(val){
			q1Res = val;
			lonlat2xyz(v0).then(function(val){
				v0Res = val;
				lonlat2xyz(v1).then(function(val){
					v1Res = val;
					quaternion(v0Res, v1Res).then(function(val){
						q2Res = val;

						quatMultiply( q1Res, q2Res).then(function(val){
							
							res(quat2euler(val));
						});
					});
				});
			});
		});	
	});	
}


/**************end of math functions**********************/