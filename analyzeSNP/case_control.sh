#!/bin/bash

sed -i "s/$2$/Case/g" $1
sed -i "s/[^Case]$/Control/g" $1
sed -i "s/Case$/1/g" $1
sed -i "s/Control$/2/g" $1