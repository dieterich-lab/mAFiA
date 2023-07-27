#!/usr/bin/env python
#
# RODAN
# v1.0
# (c) 2020,2021,2022 Don Neumann
#

import sys
import numpy as np
import torch
import torch.nn as nn

rna_default = [[-1, 256, 0, 3, 1, 1, 0], [-1, 256, 1, 10, 1, 1, 1], [-1, 256, 1, 10, 10, 1, 1], [-1, 320, 1, 10, 1, 1, 1], [-1, 384, 1, 15, 1, 1, 1], [-1, 448, 1, 20, 1, 1, 1], [-1, 512, 1, 25, 1, 1, 1], [-1, 512, 1, 30, 1, 1, 1], [-1, 512, 1, 35, 1, 1, 1], [-1, 512, 1, 40, 1, 1, 1], [-1, 512, 1, 45, 1, 1, 1], [-1, 512, 1, 50, 1, 1, 1], [-1, 768, 1, 55, 1, 1, 1], [-1, 768, 1, 60, 1, 1, 1], [-1, 768, 1, 65, 1, 1, 1], [-1, 768, 1, 70, 1, 1, 1], [-1, 768, 1, 75, 1, 1, 1], [-1, 768, 1, 80, 1, 1, 1], [-1, 768, 1, 85, 1, 1, 1], [-1, 768, 1, 90, 1, 1, 1], [-1, 768, 1, 95, 1, 1, 1], [-1, 768, 1, 100, 1, 1, 1]]
dna_default = [[-1, 320, 0, 3, 1, 1, 0], [-1, 320, 1, 3, 3, 1, 1], [-1, 384, 1, 6, 1, 1, 1], [-1, 448, 1, 9, 1, 1, 1], [-1, 512, 1, 12, 1, 1, 1], [-1, 576, 1, 15, 1, 1, 1], [-1, 640, 1, 18, 1, 1, 1], [-1, 704, 1, 21, 1, 1, 1], [-1, 768, 1, 24, 1, 1, 1], [-1, 832, 1, 27, 1, 1, 1], [-1, 896, 1, 30, 1, 1, 1], [-1, 960, 1, 33, 1, 1, 1]]

class Objectview(object):
    def __init__(self, d):
        self.__dict__ = d
        self.orig = d

def strided_app(a, L, S ):  # Window len = L, Stride len/stepsize = S
    nrows = ((a.size-L)//S)+1
    n = a.strides[0]
    return np.lib.stride_tricks.as_strided(a, shape=(nrows,L), strides=(S*n,n))


class Swish(nn.Module):
    def __init__(self, inplace=True):
        super(Swish, self).__init__()
        self.inplace = inplace

    def forward(self, x):
        return x.mul_(x.sigmoid()) if self.inplace else x.mul(x.sigmoid())


class Mish(nn.Module):
    def __init__(self):
        super().__init__()

    def forward(self, x):
        return x * (torch.tanh(torch.nn.functional.softplus(x)))


class Squeeze_Excite(torch.nn.Module):
    def __init__(self, in_channels=512, size=1, reduction="/16", activation=torch.nn.GELU):
        super(Squeeze_Excite, self).__init__()
        self.in_channels = in_channels
        self.avg = torch.nn.AdaptiveAvgPool1d(1)
        if type(reduction) == str:
            self.reductionsize = self.in_channels // int(reduction[1:])
        else:
            self.reductionsize = reduction
        self.fc1 = nn.Linear(self.in_channels, self.reductionsize)
        self.activation = activation()  # was nn.ReLU(inplace=True)
        self.fc2 = nn.Linear(self.reductionsize, self.in_channels)
        self.sigmoid = nn.Sigmoid()

    def forward(self, x):
        input = x
        x = self.avg(x)
        x = x.permute(0, 2, 1)
        x = self.activation(self.fc1(x))
        x = self.sigmoid(self.fc2(x))
        return input * x.permute(0, 2, 1)


class Convblock(torch.nn.Module):
    def __init__(self, in_channels, out_channels, kernel_size, stride=1, padding=0, dilation=1, groups=1, bias=False,
                 seperable=True, expansion=True, batchnorm=True, dropout=0.1, activation=None, sqex=True,
                 squeeze=32, sqex_activation=None, residual=True):
        # no bias?
        super(Convblock, self).__init__()
        self.seperable = seperable
        self.batchnorm = batchnorm
        self.dropout = dropout
        if activation is not None:
            self.activation = activation
        else:
            self.activation = torch.nn.GELU
        if sqex_activation is not None:
            self.sqex_activation = sqex_activation
        else:
            self.sqex_activation = torch.nn.GELU
        self.squeeze = squeeze
        self.stride = stride
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.residual = residual
        self.doexpansion = expansion
        # fix self.squeeze
        dwchannels = in_channels
        if seperable:
            if self.doexpansion and self.in_channels != self.out_channels:
                self.expansion = torch.nn.Conv1d(in_channels, out_channels, kernel_size=1, stride=1, groups=1,
                                                 bias=False)
                self.expansion_norm = torch.nn.BatchNorm1d(out_channels)
                self.expansion_act = self.activation()
                dwchannels = out_channels

            self.depthwise = torch.nn.Conv1d(dwchannels, out_channels, kernel_size=kernel_size, stride=stride,
                                             padding=padding, dilation=dilation, bias=bias,
                                             groups=out_channels // groups)
            if self.batchnorm:
                self.bn1 = torch.nn.BatchNorm1d(out_channels)
            self.act1 = self.activation()
            if self.squeeze:
                self.sqex = Squeeze_Excite(in_channels=out_channels, reduction=self.squeeze, activation=self.sqex_activation)
            self.pointwise = torch.nn.Conv1d(out_channels, out_channels, kernel_size=1, dilation=dilation, bias=bias,
                                             padding=0)
            if self.batchnorm:
                self.bn2 = torch.nn.BatchNorm1d(out_channels)
            self.act2 = self.activation()
            if self.dropout:
                self.drop = torch.nn.Dropout(self.dropout)
        else:
            self.conv = torch.nn.Conv1d(in_channels, out_channels, kernel_size=kernel_size, stride=stride,
                                        padding=padding, dilation=dilation, bias=bias)
            if self.batchnorm:
                self.bn1 = torch.nn.BatchNorm1d(out_channels)
            self.act1 = self.activation()
            if self.squeeze:
                self.sqex = Squeeze_Excite(in_channels=out_channels, reduction=self.squeeze, activation=self.sqex_activation)
            if self.dropout:
                self.drop = torch.nn.Dropout(self.dropout)
        if self.residual and self.stride == 1:
            self.rezero = nn.Parameter(torch.Tensor([0]), requires_grad=True)

    def forward(self, x):
        orig = x

        if self.seperable:
            if self.in_channels != self.out_channels and self.doexpansion:
                x = self.expansion(x)
                x = self.expansion_norm(x)
                x = self.expansion_act(x)
            x = self.depthwise(x)
            if self.batchnorm: x = self.bn1(x)
            x = self.act1(x)
            if self.squeeze:
                x = self.sqex(x)
            x = self.pointwise(x)
            if self.batchnorm: x = self.bn2(x)
            x = self.act2(x)
            if self.dropout: x = self.drop(x)
        else:
            x = self.conv(x)
            if self.batchnorm: x = self.bn1(x)
            x = self.act1(x)
            if self.dropout: x = self.drop(x)

        if self.residual and self.stride == 1 and self.in_channels == self.out_channels and x.shape[2] == orig.shape[2]:
            return orig + self.rezero * x  # rezero
            # return orig + x # normal residual
        else:
            return x


def Activation_Function(activation):
    if activation == "mish":
        return Mish
    elif activation == "swish":
        return Swish
    elif activation == "relu":
        return torch.nn.ReLU
    elif activation == "gelu":
        return torch.nn.GELU
    else:
        print("Unknown activation type:", activation)
        sys.exit(1)


class Rodan(nn.Module):
    def __init__(self, config=None, arch=None, debug=False):
        super().__init__()
        if debug: print("Initializing network")

        self.seqlen = config.seqlen
        self.vocab = config.vocab
        self.bn = nn.BatchNorm1d
        arch = rna_default if arch is None else arch

        activation = Activation_Function(config.activation.lower())
        sqex_activation = Activation_Function(config.sqex_activation.lower())

        self.convlayers = nn.Sequential()
        in_channels = 1
        convsize = self.seqlen

        for i, layer in enumerate(arch):
            paddingarg, out_channels, seperable, kernel, stride, sqex, dodropout = layer

            expansion = True

            if dodropout:
                dropout = config.dropout
            else:
                dropout = 0
            if sqex:
                squeeze = config.sqex_reduction
            else:
                squeeze = 0

            if paddingarg == -1:
                padding = kernel // 2
            else:
                padding = paddingarg
            if i == 0: expansion = False

            convsize = (convsize + (padding * 2) - (kernel - stride)) // stride
            if debug:
                print("padding:", padding, "seperable:", seperable, "ch", out_channels, "k:", kernel, "s:", stride,
                      "sqex:", sqex, "drop:", dropout, "expansion:", expansion)
                print("convsize:", convsize)
            self.convlayers.add_module("conv" + str(i),
                                       Convblock(in_channels, out_channels, kernel, stride=stride, padding=padding,
                                                 seperable=seperable, activation=activation, expansion=expansion,
                                                 dropout=dropout, squeeze=squeeze, sqex_activation=sqex_activation,
                                                 residual=True))
            in_channels = out_channels
            self.final_size = out_channels

        self.final = nn.Linear(self.final_size, len(self.vocab))
        if debug: print("Finished init network")

    def forward(self, x):
        x = self.convlayers(x)
        x = x.permute(0, 2, 1)
        x = self.final(x)

        x = torch.nn.functional.log_softmax(x, 2)
        return x.permute(1, 0, 2)