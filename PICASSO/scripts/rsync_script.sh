#!/bin/bash
mount | grep 'backup'
if (( $? )); then
   echo backup not mounted
   mount /media/backup
   if (( $? )); then
      echo error mounting backup drive
      exit;
   fi
fi

mount | grep 'master'
if (( $? )); then
   echo master not mounted
   mount /media/master
   if (( $? )); then
      echo error mounting master drive
      exit;
   fi
fi

rsync -av /media/master/ /media/backup/
