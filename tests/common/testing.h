/*
 This Source Code Form is subject to the terms of the Mozilla Public
 License, v. 2.0. If a copy of the MPL was not distributed with this
 file, You can obtain one at http://mozilla.org/MPL/2.0/. 
*/

// Part of this file has been adapted from the Eigen library.


#ifndef _PATATE_TESTING_H_
#define _PATATE_TESTING_H_

#include <iostream>
#include <vector>
#include <cerrno>
#include <cstdlib>
#include <sstream>

#define DEFAULT_REPEAT 10

#define PATATE_PP_MAKE_STRING2(S) #S
#define PATATE_PP_MAKE_STRING(S) PATATE_PP_MAKE_STRING2(S)

static std::vector<std::string> g_test_stack;
static int g_repeat;
static unsigned int g_seed;
static bool g_has_set_repeat, g_has_set_seed;

void verify_impl(bool condition, const char *testname, const char *file, int line, const char *condition_as_string)
{
  if (!condition)
  {
    std::cerr << "Test " << testname << " failed in " << file << " (" << line << ")"
      << std::endl << "    " << condition_as_string << std::endl;
    abort();
  }
}

#define VERIFY(a) ::verify_impl(a, g_test_stack.back().c_str(), __FILE__, __LINE__, PATATE_PP_MAKE_STRING(a))

#define CALL_SUBTEST(FUNC) do { \
    g_test_stack.push_back(PATATE_PP_MAKE_STRING(FUNC)); \
    FUNC; \
    g_test_stack.pop_back(); \
  } while (0)

inline void set_repeat_from_string(const char *str)
{
  errno = 0;
  g_repeat = int(strtoul(str, 0, 10));
  if(errno || g_repeat <= 0)
  {
    std::cout << "Invalid repeat value " << str << std::endl;
    exit(EXIT_FAILURE);
  }
  g_has_set_repeat = true;
}

inline void set_seed_from_string(const char *str)
{
  errno = 0;
  g_seed = int(strtoul(str, 0, 10));
  if(errno || g_seed == 0)
  {
    std::cout << "Invalid seed value " << str << std::endl;
    exit(EXIT_FAILURE);
  }
  g_has_set_seed = true;
}

static bool init_testing(int argc, char *argv[])
{
  g_has_set_repeat = false;
  g_has_set_seed   = false;
  bool need_help   = false;

  for(int i = 1; i < argc; i++)
  {
    if(argv[i][0] == 'r')
    {
      if(g_has_set_repeat)
      {
        std::cout << "Argument " << argv[i] << " conflicting with a former argument" << std::endl;
        return 1;
      }
      set_repeat_from_string(argv[i]+1);
    }
    else if(argv[i][0] == 's')
    {
      if(g_has_set_seed)
      {
        std::cout << "Argument " << argv[i] << " conflicting with a former argument" << std::endl;
        return false;
      }
        set_seed_from_string(argv[i]+1);
    }
    else
    {
      need_help = true;
    }
  }

  if(need_help)
  {
    std::cout << "This test application takes the following optional arguments:" << std::endl;
    std::cout << "  rN     Repeat each test N times (default: " << DEFAULT_REPEAT << ")" << std::endl;
    std::cout << "  sN     Use N as seed for random numbers (default: based on current time)" << std::endl;
    std::cout << std::endl;
    std::cout << "If defined, the environment variables EIGEN_REPEAT and EIGEN_SEED" << std::endl;
    std::cout << "will be used as default values for these parameters." << std::endl;
    return false;
  }

  char *env_EIGEN_REPEAT = getenv("EIGEN_REPEAT");
  if(!g_has_set_repeat && env_EIGEN_REPEAT)
    set_repeat_from_string(env_EIGEN_REPEAT);
  char *env_EIGEN_SEED = getenv("EIGEN_SEED");
  if(!g_has_set_seed && env_EIGEN_SEED)
    set_seed_from_string(env_EIGEN_SEED);

  if(!g_has_set_seed) g_seed = (unsigned int) time(NULL);
  if(!g_has_set_repeat) g_repeat = DEFAULT_REPEAT;

  std::cout << "Initializing random number generator with seed " << g_seed << std::endl;
  std::stringstream ss;
  ss << "Seed: " << g_seed;
  g_test_stack.push_back(ss.str());
  srand(g_seed);
  std::cout << "Repeating each test " << g_repeat << " times" << std::endl;
  
  return true;
}

#endif // _PATATE_TESTING_H_
