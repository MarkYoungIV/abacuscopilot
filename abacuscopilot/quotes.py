"""Quotations for the farewell screen.

Each quote is a (chinese, english) pair.  ``random_quote()`` picks one;
``farewell(console)`` prints the quote and the goodbye line.
"""

from __future__ import annotations

import random

_QUOTES: list[tuple[str, str]] = [
    (
        "所有的困难都是暂时的。",
        "All difficulties are temporary.",
    ),
    (
        "当你的才华还撑不起你的野心时，那你就应该静下心来学习。",
        "When your talent cannot yet support your ambition, you should calm down and learn.",
    ),
    (
        "最曲折的路有时最简捷。",
        "The most winding road is sometimes the shortest.",
    ),
    (
        "逆境展示奇才，顺境隐没英才。",
        "Adversity reveals genius; prosperity conceals it.",
    ),
    (
        "人生如逆旅，我亦是行人。",
        "Life is like an inn we pass through; I too am but a traveler.",
    ),
    (
        "山桃红花满上头，蜀江春水拍山流。",
        "Red mountain-peach blossoms crown the heights, as the spring waters of the Shu River lap against the hills.",
    ),
    (
        "他年我若为青帝，报与桃花一处开。",
        "If one day I become the Lord of Spring, I will let the peach blossoms all bloom together.",
    ),
    (
        "多算胜，少算不胜，而况于无算乎。",
        "With much calculation one wins; with little one loses — how much more so with none at all.",
    ),
    (
        "星星之火，可以燎原。",
        "A single spark can start a prairie fire.",
    ),
    (
        "一切反动派都是纸老虎。",
        "All reactionaries are paper tigers.",
    ),
    (
        "为人民服务。",
        "Serve the people.",
    ),
    (
        "好好学习，天天向上。",
        "Study well and make progress every day.",
    ),
    (
        "世界是你们的，也是我们的，但是归根结底是你们的。",
        "The world is yours, as well as ours, but in the last analysis, it is yours.",
    ),
    (
        "实践是检验真理的唯一标准。",
        "Practice is the sole criterion for testing truth.",
    ),
    (
        "没有调查，就没有发言权。",
        "No investigation, no right to speak.",
    ),
    (
        "自己动手，丰衣足食。",
        "Do it yourself, and you'll have ample food and clothing.",
    ),
    (
        "枪杆子里面出政权。",
        "Political power grows out of the barrel of a gun.",
    ),
    (
        "人不犯我，我不犯人；人若犯我，我必犯人。",
        "We will not attack unless we are attacked; if we are attacked, we will certainly counter-attack.",
    ),
    (
        "一万年太久，只争朝夕。",
        "Ten thousand years are too long — seize the day, seize the hour.",
    ),
    (
        "下定决心，不怕牺牲，排除万难，去争取胜利。",
        "Be resolute, fear no sacrifice, and surmount every difficulty to win victory.",
    ),
    (
        "世上无难事，只要肯登攀。",
        "Nothing is hard in this world if you dare to scale the heights.",
    ),
    (
        "虚心使人进步，骄傲使人落后。",
        "Modesty helps one go forward, whereas conceit makes one lag behind.",
    ),
    (
        "多少事，从来急；天地转，光阴迫。",
        "So many deeds cry out to be done, and always urgently; the world rolls on, time presses.",
    ),
    (
        "雄关漫道真如铁，而今迈步从头越。",
        "The strong pass is a wall of iron, yet with firm strides we are crossing its summit.",
    ),
    (
        "东方不亮西方亮，黑了南方有北方。",
        "When it is dark in the east, it is light in the west; when things are dark in the south there is still light in the north.",
    ),
    (
        "军民团结如一人，试看天下谁能敌。",
        "If the army and the people are united as one, who in the world can match them?",
    ),
    (
        "自力更生，艰苦奋斗。",
        "Self-reliance and hard struggle.",
    ),
    (
        "团结就是力量。",
        "Unity is strength.",
    ),
    (
        "真理有时候在少数人手里。",
        "Truth is sometimes in the hands of a minority.",
    ),
    (
        "不怕慢，就怕站。",
        "Don't be afraid of going slowly; be afraid of standing still.",
    ),
    (
        "情况是在不断的变化，要使自己的思想适应新的情况，就得学习。",
        "Conditions are constantly changing; to adapt your thinking to new conditions, you must study.",
    ),
    (
        "把别人的经验变成自己的，他的本事就大了。",
        "He who makes others' experience his own is truly capable.",
    ),
    (
        "道路是曲折的，前途是光明的。",
        "The road is tortuous, but the future is bright.",
    ),
    (
        "革命不是请客吃饭。",
        "Revolution is not a dinner party.",
    ),
    (
        "全心全意为人民服务。",
        "Serve the people wholeheartedly.",
    ),
    (
        "一切从实际出发。",
        "Proceed from actual facts in everything.",
    ),
    (
        "没有正确的政治观点，就等于没有灵魂。",
        "Without a correct political viewpoint, one is without a soul.",
    ),
    (
        "不但要团结和自己意见相同的人，而且要善于团结那些和自己意见不同的人。",
        "Not only unite with those who share your views, but also be good at uniting with those who disagree with you.",
    ),
    (
        "在战略上要藐视敌人，在战术上要重视敌人。",
        "Strategically, despise the enemy; tactically, take the enemy seriously.",
    ),
    (
        "一切从人民的利益出发。",
        "Proceed in all things from the interests of the people.",
    ),
    (
        "不懂就是不懂，不要装懂。",
        "If you don't understand, you don't understand — don't pretend to understand.",
    ),
    (
        "看事情要看本质。",
        "Look at the essence of things.",
    ),
    (
        "人民，只有人民，才是创造世界历史的动力。",
        "The people, and the people alone, are the motive force in the making of world history.",
    ),
    (
        "我们应该谦虚谨慎，戒骄戒躁。",
        "We should be modest and prudent, guard against arrogance and rashness.",
    ),
    (
        "任何困难都难不倒英雄的中国人民。",
        "No difficulty can daunt the heroic Chinese people.",
    ),
    (
        "天若有情天亦老，人间正道是沧桑。",
        "If heaven has feelings, heaven too will grow old; the right way on earth is through the changes of the world.",
    ),
    (
        "宜将剩勇追穷寇，不可沽名学霸王。",
        "With power and to spare we must pursue the fleeing foe, and not ape Xiang Yu the conqueror seeking idle fame.",
    ),
    (
        "春风杨柳万千条，六亿神州尽舜尧。",
        "Spring winds move willow wands, in countless tufts; six hundred million in this land all equal Shun and Yao.",
    ),
    (
        "冷眼向洋看世界，热风吹雨洒江天。",
        "With a cold eye I survey the world beyond the seas; a hot wind sprinkles rain on the river and the sky.",
    ),
    (
        "坐地日行八万里，巡天遥看一千河。",
        "Sitting on earth, I travel eighty thousand li a day; surveying the sky, I see a thousand Milky Ways from afar.",
    ),
    (
        "不管风吹浪打，胜似闲庭信步。",
        "Let the wind blow and the waves beat — better far than idly strolling in a courtyard.",
    ),
    (
        "不到长城非好汉。",
        "He who has not been to the Great Wall is not a true man.",
    ),
    (
        "数风流人物，还看今朝。",
        "For truly great men, look to this very age.",
    ),
    (
        "问苍茫大地，谁主沉浮？",
        "I ask, on this vast boundless land, who decides the rise and fall?",
    ),
    (
        "恰同学少年，风华正茂。",
        "Young we were, schoolmates, at life's full flowering.",
    ),
    (
        "书生意气，挥斥方遒。",
        "Filled with student enthusiasm, boldly we cast all restraints aside.",
    ),
    (
        "红军不怕远征难，万水千山只等闲。",
        "The Red Army fears not the trials of the Long March, holding light ten thousand crags and torrents.",
    ),
    (
        "更喜岷山千里雪，三军过后尽开颜。",
        "Glad are we that the Min Mountains are capped with thousand-li snow; after the three armies cross, every face beams with joy.",
    ),
    (
        "看万山红遍，层林尽染。",
        "I see hills on hills all in red, and wood on wood in a deep dye.",
    ),
    (
        "怅寥廓，问苍茫大地，谁主沉浮？",
        "Saddened by the vastness, I ask the great earth and the boundless blue: who decides the rise and fall?",
    ),
    (
        "为有牺牲多壮志，敢教日月换新天。",
        "Bitter sacrifice strengthens bold resolve, which dares to make sun and moon shine in new skies.",
    ),
    (
        "独有英雄驱虎豹，更无豪杰怕熊罴。",
        "Only heroes can drive away tigers and leopards; no brave men fear bears.",
    ),
    (
        "四海翻腾云水怒，五洲震荡风雷激。",
        "The four seas are churning with angry clouds and waters; the five continents are shaking with storm and thunder.",
    ),
    (
        "金猴奋起千钧棒，玉宇澄清万里埃。",
        "The Golden Monkey raises his massive cudgel, and the jade-like heavens are cleared of dust for ten thousand li.",
    ),
    (
        "梅花欢喜漫天雪，冻死苍蝇未足奇。",
        "The plum blossom welcomes the whirling snow; it is no wonder that flies freeze to death.",
    ),
    (
        "凡是敌人反对的，我们就要拥护；凡是敌人拥护的，我们就要反对。",
        "We should support whatever the enemy opposes and oppose whatever the enemy supports.",
    ),
    (
        "群众是真正的英雄。",
        "The masses are the real heroes.",
    ),
    (
        "知识的问题是一个科学问题，来不得半点虚伪和骄傲。",
        "Knowledge is a matter of science, and no dishonesty or conceit whatsoever is permissible.",
    ),
    (
        "感觉到的东西，我们不能立刻理解它，只有理解了的东西才更深刻地感觉它。",
        "What we sense we cannot immediately understand; only what we understand can we more deeply sense.",
    ),
    (
        "矛盾存在于一切事物的发展过程中。",
        "Contradiction exists in the process of development of all things.",
    ),
    (
        "外因是变化的条件，内因是变化的根据。",
        "External causes are the condition of change; internal causes are the basis of change.",
    ),
    (
        "物质可以变成精神，精神可以变成物质。",
        "Matter can be transformed into spirit, and spirit into matter.",
    ),
    (
        "世界上怕就怕“认真”二字。",
        "What the world fears most is the word 'seriousness'.",
    ),
    (
        "人类的历史，就是一个不断地从必然王国向自由王国发展的历史。",
        "The history of mankind is one of continuous development from the realm of necessity to the realm of freedom.",
    ),
    (
        "青年是整个社会力量中最积极最有生气的力量。",
        "Young people are the most active and vital force in society.",
    ),
    (
        "你们青年人朝气蓬勃，正在兴旺时期，好像早晨八九点钟的太阳。",
        "You young people, full of vigor and vitality, are in the bloom of life, like the sun at eight or nine in the morning.",
    ),
    (
        "没有文化的军队是愚蠢的军队。",
        "An army without culture is a stupid army.",
    ),
    (
        "战争的目的不是别的，就是保存自己，消灭敌人。",
        "The aim of war is none other than to preserve oneself and destroy the enemy.",
    ),
    (
        "读书是学习，使用也是学习，而且是更重要的学习。",
        "Reading is learning, but applying is also learning, and indeed the more important kind.",
    ),
    (
        "凡事预则立，不预则废。",
        "Preparedness ensures success; unpreparedness spells failure.",
    ),
    (
        "错误常常是正确的先导。",
        "Error is often the precursor of what is correct.",
    ),
]


def random_quote() -> tuple[str, str]:
    """Return a random (chinese, english) quote pair."""
    return random.choice(_QUOTES)


def farewell(console) -> None:
    """Print a random Mao quote followed by the goodbye line."""
    cn, en = random_quote()
    console.print()
    console.print(f"  [dim italic]{cn}[/dim italic]")
    console.print(f"  [dim italic]{en}[/dim italic]")
    console.print()
    console.print("[yellow]Goodbye! 再见！[/yellow]")
