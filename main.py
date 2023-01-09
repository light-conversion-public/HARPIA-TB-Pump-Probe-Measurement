#!/usr/bin/env python
# -*- coding: utf-8 -*-
#==========================================================================
# 
#--------------------------------------------------------------------------
# Copyright (c) 2022 Light Conversion, UAB
# All rights reserved.
# www.lightcon.com
#==========================================================================
if __name__ == '__main__':
    import os
    import sys    
        
    sys.path.append(os.path.dirname(os.path.realpath(sys.argv[0])))
    os.chdir(os.path.dirname(os.path.realpath(sys.argv[0])))
    
    import lclauncher
        
    from package import app