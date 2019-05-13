/* -*- Mode: C++; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*-
  * ***** BEGIN LICENSE BLOCK *****
  * Version: MPL 1.1/GPL 2.0/LGPL 2.1
  *
  * The contents of this file are subject to the Mozilla Public License Version
  * 1.1 (the "License"); you may not use this file except in compliance with
  * the License. You may obtain a copy of the License at
  * http://www.mozilla.org/MPL/
  *
  * Software distributed under the License is distributed on an "AS IS" basis,
  * WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
  * for the specific language governing rights and limitations under the
  * License.
  *
  * The Original Code is Mozilla code.
  *
  * The Initial Developer of the Original Code is
  * Mozilla Foundation.
  * Portions created by the Initial Developer are Copyright (C) 2010
  * the Initial Developer. All Rights Reserved.
  *
  * Contributor(s):
  *   Taras Glek <tglek@mozilla.com>
  *
  * Alternatively, the contents of this file may be used under the terms of
  * either the GNU General Public License Version 2 or later (the "GPL"), or
  * the GNU Lesser General Public License Version 2.1 or later (the "LGPL"),
  * in which case the provisions of the GPL or the LGPL are applicable instead
  * of those above. If you wish to allow use of your version of this file only
  * under the terms of either the GPL or the LGPL, and not to allow others to
  * use your version of this file under the terms of the MPL, indicate your
  * decision by deleting the provisions above and replace them with the notice
  * and other provisions required by the GPL or the LGPL. If you do not delete
  * the provisions above, a recipient may use your version of this file under
  * the terms of any one of the MPL, the GPL or the LGPL.
  *
  * ***** END LICENSE BLOCK ***** */
 
 #include <fcntl.h>
 #include <unistd.h>
 #include <sys/types.h>
 #include <sys/stat.h>
 
 // 20150108 RJVB: created from the OSX-specific code from Mozilla's mozilla::fallocation() function
 // of which the licensing information is copied above.
 
 #ifdef cplusplus
 extern "C" {
 #endif
 
#ifdef __APPLE__
 static int posix_fallocate(int fd, off_t offset, off_t len)
 {
     fstore_t store = {F_ALLOCATECONTIG, F_PEOFPOSMODE, offset, len};
     // Try to get a continous chunk of disk space
     int ret = fcntl(fd, F_PREALLOCATE, &store);
     if(-1 == ret){
         // OK, perhaps we are too fragmented, allocate non-continuous
         store.fst_flags = F_ALLOCATEALL;
         ret = fcntl(fd, F_PREALLOCATE, &store);
         if (-1 == ret)
             return false;
     }
     return 0 == ftruncate(fd, len);
 }
#endif
 
 #ifdef cplusplus
 }
 #endif
