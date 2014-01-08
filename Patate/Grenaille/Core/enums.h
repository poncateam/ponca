/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/

#ifndef _GRENAILLE_ENUMS_
#define _GRENAILLE_ENUMS_

namespace Grenaille
{
	/*! Enum corresponding to the state of a fitting method (and what the finalize function can return) */
	enum FITRESULT
	{
		STABLE = 0, /*!< The fitting is stable an ready to use (and having more than 6 neighbours)*/
		UNSTABLE = 1, /*!< The fitting is ready to use but it can be unstable (and having between 3 and 6 neighbors)*/
		UNDEFINED = 2, /*!< The fitting is undefined, you can't use it for valid results (and having less than 3 neighbors)*/
		NBMAX /*!< Nb enums */
	};

} //namespace Grenaille

#endif
