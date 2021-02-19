#!/usr/bin/env python

import click
import SimpleITK as sitk
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import scipy.special as special


@click.command()
@click.option("-i", "--input", required=True, help="Input file name")
@click.option("-o", "--output", required=True, help="Output file name")
@click.option("-m", "--m", required=True, default=500, help="Size of the reconstruction image matrix mxm")
@click.option("-k", "--k",required=True, default=500, help="Size of the region used for offset correction kxk")

def BPF(input, output, m, k):

    def BPFfilter(N, dx):
        x = y = np.arange(N) * dx - (N / 2 - 1) *dx #filter sampled such that x=y=0 is a sample otherwise artifacts appear
        X , Y = np.meshgrid (x, y)
        k = np.zeros(X.shape)
        R = np.sqrt(X**2 + Y**2)
        not_zero = R!=0
        var_k = np.pi * R[not_zero] / dx
        k[not_zero] = 1./(4 * np.pi**2 * R[not_zero]**3) * (var_k**2 * special.jn(1,var_k) - np.pi * var_k/2 * (special.jn(1,var_k) * special.struve(0,var_k) - special.jn(0,var_k) * special.struve(1,var_k)))
        k[R==0] = np.pi / (12 * dx**3)
        return k

    corrSize = int(k)
    reconSize = int(m)
    imgReader = sitk.ImageFileReader()
    imgReader.SetFileName(input)
    sino = imgReader.Execute()
    spacing = sino.GetSpacing()
    offset = sino.GetOrigin()
    sino = sitk.GetArrayFromImage(sino)
    bpSize = sino.shape[1]
    nb_proj = sino.shape[0]
    dphi = np.pi / nb_proj
    proj_axis = 0
    dx = spacing[-1]

    #Compute backprojection
    bp = np.sum(sino, axis = proj_axis) * dphi
    sino = []

    #Compute filter on very large matrix
    k = BPFfilter(corrSize, dx)
    a = int((corrSize - bpSize) / 2) #to resize filter to BP size

    #Calculations for offset correction
    xbp = np.arange(bpSize) * dx - (bpSize - 1) / 2 * dx
    Xbp, Ybp = np.meshgrid (xbp, xbp)
    xc = np.arange(corrSize) * dx - (corrSize - 1) / 2 * dx
    Xc, Yc = np.meshgrid (xc, xc)
    idx = np.in1d(Xc, Xbp) * np.in1d(Yc, Ybp)
    idx = np.reshape(idx, Xc.shape)
    xcind = np.arange(corrSize)
    Xcind, Ycind = np.meshgrid(xcind, xcind)
    SecondIntegral = (dx**2 * np.sum(k[corrSize - 1 - Xcind[~idx], corrSize - 1 - Ycind[~idx]] / (np.sqrt(Xc[~idx]**2 + Yc[~idx]**2))))
    b = int((bpSize - reconSize) / 2) #resize bp to recon size

    reconstruction = np.zeros(bp.shape)
    for z in np.arange(bp.shape[1]):
        reconstruction[:,z,:] = signal.fftconvolve(bp[:,z,:], k[a:-a,a:-a], mode='same') * dx**2
        delta = (dx**2 * np.sum(reconstruction[b:-b,z,b:-b])) * SecondIntegral
        reconstruction[:,z,:] += delta

    offset_recon = -(reconSize - 1) / 2 * dx
    reconstruction = reconstruction[b:-b,:,b:-b]
    reconimg = sitk.GetImageFromArray(reconstruction)
    reconimg.SetSpacing(spacing[:-1])
    reconimg.SetOrigin([offset_recon, offset[1],offset_recon])
    writer = sitk.ImageFileWriter()
    writer.SetFileName(output)
    writer.Execute(reconimg)

if __name__ == '__main__':
    BPF()
