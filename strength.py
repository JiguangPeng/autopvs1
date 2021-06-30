#!/usr/bin/env python
# -*- coding:utf-8 -*-
# author: Jiguang Peng
# created: 2019/6/27 16:48

from enum import Enum


class Strength(Enum):
    """
    ACMG pathogenic strength
    """
    Unset = -1
    Unmet = 0
    Supporting = 1
    Moderate = 2
    Strong = 3
    VeryStrong = 4

    def upgrade(self, level):
        if self.value + level <= 4:
            return Strength(self.value + level)
        else:
            return Strength(4)

    def downgrade(self, level):
        if self.value - level >= 0:
            return Strength(self.value - level)
        else:
            return Strength(0)
