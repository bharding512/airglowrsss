//=============================================================================
// File name  : ksbigusbmain.h                                                   
// Description: The main header file of the USB kernel driver.			     
// Author     : Jan Soldan    - Linux                   
// Author     : Matt Longmire - Windows 	     
// Copyright (C) 2003 Jan Soldan, Matt Longmire - SBIG			     
// All rights reserved.							     
//=============================================================================
#ifndef _KSBIGUSBMAIN_H_
#define _KSBIGUSBMAIN_H_
//=============================================================================
#include <linux/config.h>
#include <linux/kernel.h>
#include <linux/sched.h>
#include <linux/signal.h>
#include <linux/errno.h>
#include <linux/poll.h>
#include <linux/init.h>
#include <linux/slab.h>
#include <linux/fcntl.h>
#include <linux/module.h>
#include <linux/spinlock.h>
#include <linux/list.h>
#include <linux/smp_lock.h>
#include <linux/devfs_fs_kernel.h>
#include <linux/usb.h>	

#include "sbigdef.h"
#include "ksbigictl.h"
	    
#define	USB_BULK_MSG_TIMEOUT  (10 * HZ) 
//=============================================================================
struct device_private_data{
 struct usb_device     *device;		  
 struct usb_interface  *interface;	  
 devfs_handle_t		devfs;	  
 unsigned char		subminor;	  
 int			open_count;	 
 struct semaphore	semaphore; 	 

 struct urb            *in_urb;		  
 unsigned char         *in_buffer;	   
 int			in_size;	  
 __u8			in_endpointAddr;  
					 
 struct urb            *out_urb;	   
 unsigned char         *out_buffer;	   
 int			out_size;	   
 __u8			out_endpointAddr;  
};
//=============================================================================
int   KSbigUsbModuleInit(void);
void  KSbigUsbModuleExit(void);

void *KSbigUsbProbeDevice(
      struct usb_device          *usbDevice,
      unsigned int                interfaceNumber,
      const struct usb_device_id *usbDeviceIdTable);

void  KSbigUsbDisconnectDevice(
      struct usb_device          *usbDevice,
      void                       *drvContext);

void  KSbigUsbDeviceDelete(struct device_private_data *dpd);
 
void  KSbigUsbGetDeviceInfo(
      struct usb_device          *usbDevice,
      unsigned int                interfaceNumber,
      const struct usb_device_id *usbDeviceIdTable);

int   KSbigUsbOpen(
      struct inode               *inode,
      struct file                *filp);

int   KSbigUsbRelease(
      struct inode               *inode,
      struct file                *filp);

int   KSbigUsbIoctl(
      struct inode	         *inode,
      struct file		 *filp,
      unsigned int	          cmd,
      unsigned long	          arg);

void  KSbigUsbOutBulkCallback(
      struct urb                 *urb);

void  KSbigUsbInBulkCallback(
      struct urb                 *urb);

int   KSbigUsbTestCommand(void);

int   KSbigUsbGetDriverInfo(
      GetDriverInfoResults0      *results);

int   KSbigUsbSendMicroBlock(     
      struct device_private_data *dpd,
      LinuxMicroblock            *arg);

int   KSbigUsbGetMicroBlock(
      struct device_private_data *dpd, 
      LinuxMicroblock            *arg);

int   KSbigUsbSubmitInUrb(
      struct device_private_data *dpd, 
      LinuxMicroblock            *arg);

int   KSbigUsbGetInUrb(
      struct device_private_data *dpd, 
      LinuxMicroblock            *arg);

int   KSbigUsbGetJiffies(
      unsigned long              *arg);

int   KSbigUsbGetHz(
      unsigned long              *arg);

//=============================================================================
#endif // _KSBIGUSBMAIN_H_




