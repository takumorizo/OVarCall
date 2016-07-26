#!/usr/bin/env python


class WarningTypeException(Exception): pass
class CliticalTypeException(Exception): pass

class InvalidPileupReference(WarningTypeException):pass
class InvalidPileupFormat(CliticalTypeException):pass
class InvalidBamFileManagement(CliticalTypeException): pass
class InvalidMutationType(WarningTypeException): pass
class InsufficientBamReadFilter(WarningTypeException): pass
class TooMuchWindowSizeAllocation(CliticalTypeException): pass
class TooDeepBamPositionException(WarningTypeException): pass