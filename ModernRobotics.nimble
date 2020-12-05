# Package

version       = "0.1.0"
author        = "nimbotics"
description   = "the Modern Robotics Library rewritten in nim"
license       = "MIT"
srcDir        = "src"


# Dependencies

requires "nim >= 1.2.0"
requires "neo >= 0.3.1"


# testing

task test, "Runs the test suite":
  exec "nim c -r tests/test1"
